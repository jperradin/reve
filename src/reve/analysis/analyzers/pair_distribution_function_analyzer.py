import numpy as np
import os
from itertools import combinations
from typing import Dict, Optional, List
from scipy.spatial import cKDTree

from .base_analyzer import BaseAnalyzer
from ...core.frame import Frame
from ...config.settings import Settings
from ...utils.geometry import calculate_pbc_distance, cartesian_to_fractional


class PairDistributionFunctionAnalyzer(BaseAnalyzer):
    """
    Computes the pair distribution function G(r) for each pair in the system
    using an efficient cKDTree-based method.
    """

    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self.r: Optional[np.ndarray] = None
        self.gr: Optional[Dict[str, np.ndarray]] = None
        self._atoms_data: Optional[Dict[str, np.ndarray]] = None
        self.gr_data: list = []
        self.frame_count: int = 0
        self.r_max: float = 10.0
        self.bins: int = 800
        self.dr: float = self.r_max / self.bins

    def _check_rmax(self, lattice: np.ndarray) -> None:
        """Ensures r_max is not larger than half the box size."""
        min_box_dim = np.min(np.linalg.norm(lattice, axis=1))
        if self.r_max > min_box_dim / 2.0:
            self.r_max = min_box_dim / 2.0
            self.bins = int(self.r_max / self.dr)
            if self._settings.verbose:
                print(
                    f"Warning: r_max for PDF has been rescaled to half the smallest box dimension: {self.r_max:.3f}"
                )

    def analyze(self, frame: Frame) -> None:
        self._atoms_data = frame.nodes_data.wrapped_positions
        self._correlation_lengths: Dict[str, float] = (
            frame.nodes_data.correlation_lengths
        )
        self._check_rmax(lattice=frame.get_lattice())
        pairs = self._get_pairs(list(self._atoms_data.keys()))
        self.calculate_gr(pairs, frame)
        if self.gr is not None:
            self.gr_data.append(self.gr)
            self.frame_count += 1

    def finalize(self) -> None:
        if self.frame_count == 0:
            return

        # Initialize a dictionary to hold the sum of G(r)
        gr_sum = {
            key: np.zeros_like(self.gr_data[0][key]) for key in self.gr_data[0].keys()
        }

        # Sum the G(r) from all frames
        for gr_frame in self.gr_data:
            for key, value in gr_frame.items():
                gr_sum[key] += value

        # Calculate the average
        self.gr = {key: value / self.frame_count for key, value in gr_sum.items()}

    def get_result(self) -> Dict[str, float]:
        return self.gr if self.gr is not None else {}

    def print_to_file(self) -> None:
        self.finalize()
        if self.gr is None or self.r is None:
            return
        output_path = os.path.join(
            self._settings.export_directory, "pair_distribution_function.dat"
        )
        with open(output_path, "w") as f:
            header = f"# r\t" + "\t".join(self.gr.keys()) + "\n"
            f.write(header)
            for i in range(len(self.r)):
                line = f"{self.r[i]:.5f}"
                for key in self.gr.keys():
                    line += f"\t{self.gr[key][i]:.5f}"
                f.write(line + "\n")

    def _get_pairs(self, species: list[str]) -> list[str]:
        pairs = [f"{s1}-{s2}" for s1, s2 in combinations(species, 2)]
        for s in species:
            pairs.append(f"{s}-{s}")
        return pairs

    def calculate_gr(self, pairs: list[str], frame: Frame) -> None:
        lattice = frame.get_lattice()
        box_volume = np.linalg.det(lattice)
        inv_lattice = np.linalg.inv(lattice)

        self.r = np.linspace(0, self.r_max, self.bins + 1)
        r_mid = (self.r[:-1] + self.r[1:]) / 2
        shell_volumes = 4.0 / 3.0 * np.pi * (self.r[1:] ** 3 - self.r[:-1] ** 3)

        self.gr = {}

        # Use fractional coordinates for robust PBC handling in cKDTree
        positions_frac = {
            s: cartesian_to_fractional(p, lattice) for s, p in self._atoms_data.items()
        }

        # A safe search radius in fractional coordinates to find all candidate pairs
        search_radius_frac = self.r_max * np.linalg.norm(inv_lattice, 2)

        for pair in pairs:
            s1, s2 = pair.split("-")
            pos1, pos2 = self._atoms_data[s1], self._atoms_data[s2]
            pos1_frac, pos2_frac = positions_frac[s1], positions_frac[s2]
            n1, n2 = len(pos1), len(pos2)

            if n1 == 0 or n2 == 0:
                continue

            tree2 = cKDTree(pos2_frac, boxsize=[1, 1, 1])
            distances = []

            if s1 == s2:
                if n1 < 2:
                    continue
                # Find pairs within the same species
                candidate_pairs = tree2.query_pairs(
                    r=search_radius_frac, output_type="ndarray"
                )
                for i, j in candidate_pairs:
                    dist = calculate_pbc_distance(pos1[i], pos1[j], lattice)
                    if dist < self.r_max:
                        distances.append(dist)
            else:
                # Find pairs between two different species
                tree1 = cKDTree(pos1_frac, boxsize=[1, 1, 1])
                candidate_indices = tree1.query_ball_tree(tree2, r=search_radius_frac)
                for i, indices_j in enumerate(candidate_indices):
                    for j in indices_j:
                        dist = calculate_pbc_distance(pos1[i], pos2[j], lattice)
                        if dist < self.r_max:
                            distances.append(dist)

            # Bin the collected distances
            hist, _ = np.histogram(distances, bins=self.r)

            # --- Normalization ---
            g_r = np.zeros(self.bins)
            if s1 == s2:
                norm_factor = (n1 * (n1 - 1) / 2.0) / box_volume
            else:
                norm_factor = (n1 * n2) / box_volume

            ideal_counts = shell_volumes * norm_factor
            # Avoid division by zero for bins with no ideal counts
            non_zero_mask = ideal_counts > 0
            g_r[non_zero_mask] = hist[non_zero_mask] / ideal_counts[non_zero_mask]

            self.gr[pair] = g_r

        # Calculate the total pair distribution function with Faber-Ziman formalism
        self.gr["total"] = np.zeros_like(self.gr[pairs[0]])

        cl = self._correlation_lengths  # These are the neutron scattering lengths
        species_list = sorted(self._atoms_data.keys())
        species_counts = {s: len(self._atoms_data[s]) for s in species_list}
        num_atoms = sum(species_counts.values())

        concentrations = {s: species_counts[s] / num_atoms for s in species_list}

        # Correct normalization for Faber-Ziman
        b_avg = sum(concentrations[s] * cl[s] for s in species_list)
        b_avg_gr = b_avg**2

        for s1 in species_list:
            self.gr["total"] += (
                concentrations[s1] ** 2 * cl[s1] ** 2 * self.gr[f"{s1}-{s1}"]
            )

        for s1, s2 in combinations(species_list, 2):
            pair_key = f"{s1}-{s2}"
            self.gr["total"] += (
                2
                * concentrations[s1]
                * concentrations[s2]
                * cl[s1]
                * cl[s2]
                * self.gr[pair_key]
            )

        if b_avg_gr > 0:
            self.gr["total"] /= b_avg_gr

        # Update self.r to be the bin centers for plotting
        self.r = r_mid
