# src/reve/analysis/analyzers/neutron_structure_factor_fft_analyzer.py

import numpy as np
import os
from dataclasses import dataclass
from itertools import combinations
from tqdm import tqdm

from .base_analyzer import BaseAnalyzer
from ...core.frame import Frame
from ...config.settings import Settings

@dataclass
class QVectors:
    """Holds reciprocal spaces q-vectors and their binning informations."""
    q_norm: np.ndarray
    q_bins: np.ndarray

class NeutronStructureFactorFFTAnalyzer(BaseAnalyzer):
    """
    Computes the neutron structure factor using an efficient FFT-based method.
    """

    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self.q = None
        self.nsf = None
        self._atoms_data = None
        # For optimal FFT performance, grid size should be a power of 2
        self.grid_size = 128

    def analyze(self, frame: Frame) -> None:
        self._atoms_data = frame.get_wrapped_positions_by_element()
        self._correlation_lengths = frame.get_correlation_lengths()
        pairs = self._get_pairs(self._atoms_data.keys())
        self.calculate_s_q_fft(pairs, frame)

    def finalize(self) -> None:
        pass

    def get_result(self) -> dict[str, float]:
        return self.nsf

    def print_to_file(self) -> None:
        if self.nsf is None:
            return
        output_path = os.path.join(
            self._settings.export_directory, "neutron_structure_factor_fft.dat"
        )
        with open(output_path, "w") as f:
            header = "# q\t" + "\t".join(self.nsf.keys()) + "\n"
            f.write(header)
            for i in range(len(self.q)):
                line = f"{self.q[i]:.5f}"
                for key in self.nsf.keys():
                    line += f"\t{self.nsf[key][i]:.5f}"
                f.write(line + "\n")

    def _get_pairs(self, species: list[str]) -> list[str]:
        pairs = [f"{s1}-{s2}" for s1, s2 in combinations(species, 2)]
        for s in species:
            pairs.append(f"{s}-{s}")
        pairs.append("total")
        return pairs

    def _progress_iterator(self, iterable, desc: str, unit: str = "it"):
        """Returns a tqdm iterator if not in quiet mode, otherwise returns the original iterable."""
        if self._settings.verbose:
            return tqdm(iterable, desc=desc, colour="#00ffff", leave=False, unit=unit)
        return iterable

    def calculate_s_q_fft(self, pairs: list[str], frame: Frame) -> None:
        """
        Calculates the neutron structure factor using an FFT-based method.
        """
        grid_shape = (self.grid_size, self.grid_size, self.grid_size)
        lattice = frame.get_lattice()

        rho_q = {}
        iterator = self._progress_iterator(self._atoms_data.items(), "Processing species", "species")
        for species, positions in iterator:
            density_grid = self._create_density_grid(positions, lattice, grid_shape)
            rho_q[species] = np.fft.fftn(density_grid)

        q_vectors = self._generate_fft_q_vectors(lattice, grid_shape)
        f_unbinned = self._calculate_factors_from_rho_q(pairs, rho_q)
        structure_factor_binned = self._bin_structure_factor(f_unbinned, q_vectors)

        self.q = structure_factor_binned.pop("q")
        self.nsf = structure_factor_binned

    def _create_density_grid(self, positions: np.ndarray, lattice: np.ndarray, grid_shape: tuple) -> np.ndarray:
        """Bins atomic positions onto a 3D grid."""
        density_grid = np.zeros(grid_shape, dtype=np.float32)
        inv_lattice = np.linalg.inv(lattice)
        
        # Transform cartesian coordinates to grid indices
        frac_coords = (positions @ inv_lattice) % 1.0
        grid_indices = (frac_coords * np.array(grid_shape)).astype(int)

        np.add.at(density_grid, tuple(grid_indices.T), 1)
        return density_grid

    def _generate_fft_q_vectors(self, lattice: np.ndarray, grid_shape: tuple) -> QVectors:
        """Generates the q-vectors corresponding to the FFT grid."""
        q_axes = [2 * np.pi * np.fft.fftfreq(n, d=1.0) for n in grid_shape]
        q_grid = np.stack(np.meshgrid(*q_axes, indexing='ij'), axis=-1)
        
        # Transform from fractional reciprocal coordinates to cartesian
        recip_lattice = np.linalg.inv(lattice).T
        q_vecs_cartesian = q_grid @ recip_lattice
        q_norm = np.linalg.norm(q_vecs_cartesian, axis=-1)
        
        q_max = 10.0
        q_bins = np.linspace(0, q_max, num=100) # Define bins for averaging
        return QVectors(q_norm=q_norm, q_bins=q_bins)

    def _calculate_factors_from_rho_q(self, pairs: list[str], rho_q: dict) -> dict:
        """Calculates partial and total structure factors from Fourier-space densities."""
        f = {}
        species_list = sorted(self._atoms_data.keys())
        species_counts = {s: len(self._atoms_data[s]) for s in species_list}
        num_atoms = sum(species_counts.values())

        for pair in self._progress_iterator(pairs, "Calculating partial factors", "pair"):
            if pair == "total": continue
            s1, s2 = pair.split('-')
            
            if s1 == s2:
                f[pair] = np.abs(rho_q[s1])**2 / num_atoms
            else:
                f[pair] = (rho_q[s1] * np.conj(rho_q[s2])).real / num_atoms
        
        if "total" in pairs:
            f["total"] = self._calculate_total_factor(f, species_list, species_counts, num_atoms)
        return f

    def _calculate_total_factor(self, f_partial, species_list, species_counts, num_atoms):
        """Calculates the total structure factor from partial factors."""
        f_total = np.zeros_like(next(iter(f_partial.values())))
        cl = self._correlation_lengths

        for i, species in enumerate(species_list):
            f_total += cl[species] ** 2 * species_counts[species] / num_atoms * f_partial[f"{species}-{species}"]

        for s1, s2 in combinations(species_list, 2):
            pair_key = f"{s1}-{s2}"
            f_total += 2 * cl[s1] * cl[s2] * np.sqrt(species_counts[s1] * species_counts[s2]) / num_atoms * f_partial[pair_key]

        normalization = sum(species_counts[s] * cl[s] ** 2 for s in species_list)
        f_total *= num_atoms / normalization
        return f_total

    def _bin_structure_factor(self, f_unbinned: dict, q_vectors: QVectors) -> dict:
        """Vectorized binning of structure factors by averaging over shells in q-space."""
        q_norm_flat = q_vectors.q_norm.ravel()
        q_bins = q_vectors.q_bins
        num_bins = len(q_bins) - 1

        is_not_q_zero_mask = q_norm_flat > 1e-9
        q_norms_to_bin = q_norm_flat[is_not_q_zero_mask]

        bin_indices = np.digitize(q_norms_to_bin, q_bins) - 1
        is_in_valid_bin_mask = (bin_indices >= 0) & (bin_indices < num_bins)
        final_bin_indices = bin_indices[is_in_valid_bin_mask]

        q_midpoints = (q_bins[:-1] + q_bins[1:]) / 2
        binned_sf = {"q": q_midpoints}

        counts = np.bincount(final_bin_indices, minlength=num_bins)
        counts[counts == 0] = 1

        for pair, unbinned_values in f_unbinned.items():
            sf_to_bin = unbinned_values.ravel()[is_not_q_zero_mask]
            sf_weights = sf_to_bin[is_in_valid_bin_mask]
            sums = np.bincount(final_bin_indices, weights=sf_weights, minlength=num_bins)
            binned_sf[pair] = sums / counts
            
        return binned_sf