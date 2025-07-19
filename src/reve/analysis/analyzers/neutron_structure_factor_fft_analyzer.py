# src/reve/analysis/analyzers/neutron_structure_factor_fft_analyzer.py

import numpy as np
import os
from dataclasses import dataclass
from itertools import combinations
from tqdm import tqdm
from scipy.stats import binned_statistic

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
        self.q_max = 10.0
        self.q_bins = 100

    def analyze(self, frame: Frame) -> None:
        self._atoms_data = frame.get_wrapped_positions_by_element()
        self._correlation_lengths = frame.get_correlation_lengths()
        pairs = self._get_pairs(list(self._atoms_data.keys()))
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
            for i in range(2,len(self.q)):
                # Skipping 2 values
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
        box_size = np.diag(lattice)

        rho_q = {}
        iterator = self._progress_iterator(self._atoms_data.items(), "Processing species", "species")
        for species, positions in iterator:
            density_grid, _ = np.histogramdd(positions, bins=grid_shape, range=[(0, L) for L in box_size])
            rho_q[species] = np.fft.fftn(density_grid)

        q_vectors = self._generate_fft_q_vectors(lattice, grid_shape)
        f_unbinned = self._calculate_factors_from_rho_q(pairs, rho_q)
        
        # Shift the zero-frequency component to the center for easier radial binning
        for pair in f_unbinned:
            f_unbinned[pair] = np.fft.fftshift(f_unbinned[pair])
        q_vectors.q_norm = np.fft.fftshift(q_vectors.q_norm)

        structure_factor_binned = self._bin_structure_factor(f_unbinned, q_vectors)

        self.q = structure_factor_binned.pop("q")
        self.nsf = structure_factor_binned

    def _generate_fft_q_vectors(self, lattice: np.ndarray, grid_shape: tuple) -> QVectors:
        """Generates the q-vectors corresponding to the FFT grid."""
        box_size = np.diag(lattice)
        q_vecs = [np.fft.fftfreq(n, d=L/n) * 2 * np.pi for n, L in zip(grid_shape, box_size)]
        qx, qy, qz = np.meshgrid(q_vecs[0], q_vecs[1], q_vecs[2], indexing='ij')
        q_norm = np.sqrt(qx**2 + qy**2 + qz**2)
        
        q_bins = np.linspace(0, self.q_max, self.q_bins)
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
                if species_counts[s1] > 0:
                    f[pair] = (np.abs(rho_q[s1])**2) / species_counts[s1]
                else:
                    f[pair] = np.zeros_like(rho_q[s1])
            else:
                if species_counts[s1] > 0 and species_counts[s2] > 0:
                    f[pair] = (rho_q[s1] * np.conj(rho_q[s2])).real / np.sqrt(species_counts[s1] * species_counts[s2])
                else:
                    f[pair] = np.zeros_like(rho_q[s1])
        
        if "total" in pairs:
            f["total"] = self._calculate_total_factor(f, species_list, species_counts, num_atoms)

        return f

    def _calculate_total_factor(self, f_partial, species_list, species_counts, num_atoms):
        """Calculates the total structure factor from partial factors using the Faber-Ziman formalism."""
        f_total = np.zeros_like(next(iter(f_partial.values())))
        cl = self._correlation_lengths # These are the neutron scattering lengths
        
        concentrations = {s: species_counts[s] / num_atoms for s in species_list}
        b_avg_sq = sum(concentrations[s] * cl[s]**2 for s in species_list)

        for s1 in species_list:
            f_total += concentrations[s1]**2 * cl[s1]**2 * f_partial[f"{s1}-{s1}"]
        
        for s1, s2 in combinations(species_list, 2):
            pair_key = f"{s1}-{s2}"
            f_total += 2 * concentrations[s1] * concentrations[s2] * cl[s1] * cl[s2] * f_partial[pair_key]

        if b_avg_sq > 0:
            f_total /= b_avg_sq

        return f_total

    def _bin_structure_factor(self, f_unbinned: dict, q_vectors: QVectors) -> dict:
        """Vectorized binning of structure factors by averaging over shells in q-space."""
        q_norm_flat = q_vectors.q_norm.ravel()
        q_bins = q_vectors.q_bins
        
        q_binned = (q_bins[:-1] + q_bins[1:]) / 2
        binned_sf = {"q": q_binned}

        for pair, unbinned_values in f_unbinned.items():
            binned_values, _, _ = binned_statistic(
                q_norm_flat,
                unbinned_values.ravel(),
                statistic='mean',
                bins=q_bins
            )
            binned_sf[pair] = np.nan_to_num(binned_values) # Replace NaNs from empty bins with 0
            
        return binned_sf