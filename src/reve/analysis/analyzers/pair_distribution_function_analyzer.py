import numpy as np
import os
from tqdm import tqdm
from itertools import combinations
from numba import njit, prange
from typing import Dict, Optional, List

from .base_analyzer import BaseAnalyzer
from ...core.frame import Frame
from ...config.settings import Settings
from ...utils.geometry import calculate_pbc_distance


class PairDistributionFunctionAnalyzer(BaseAnalyzer):
    """
    Computes the pair distribution function G(r) for each pair in the system.
    """

    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self.r = None
        self.gr = None
        self._atoms_data: Optional(Dict[str, nd.ndarray]) = None
        self.gr_data = []
        self.frame_count = 0
        self.r_max = 10.0
        self.bins = 600
        self.dr = float(self.r_max / self.bins)

        print("init pair_distribution_function")

    def _check_rmax(self, lattice: np.ndarray) -> None:
        check = False
        for l in np.diag(lattice):
            if self.r_max > l / 2:
                check = True
        if check:
            self.r_max = min(np.diag(lattice)) / 2
            print(
                f"Warning: r_max has been rescaled to half the smallest box dimension: {self.r_max}"
            )

    def analyze(self, frame: Frame) -> None:
        self._atoms_data = frame.nodes_data.wrapped_positions
        self._check_rmax(lattice=frame.get_lattice())
        pairs = self._get_pairs(list(self._atoms_data.keys()))
        self.calculate_gr(pairs, frame)
        self.gr_data.append(self.gr)
        self.frame_count += 1

    def finalize(self) -> None:
        if self.frame_count == 0:
            return

        # Initialyze a dictionary to hold the sum of G(r)
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
        return self.gr

    def print_to_file(self) -> None:
        if self.gr is None:
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
        box_volume = np.linalg.det(frame.get_lattice())
        total_atoms = sum(len(pos) for pos in self._atoms_data.values())
        number_density = total_atoms / box_volume

        self.r = np.linspace(0, self.r_max, self.bins)
        gr = {}

        for pair in pairs:
            s1, s2 = pair.split("-")
            pos1 = self._atoms_data[s1]
            pos2 = self._atoms_data[s2]

            hist = self._calculate_histogram(
                pos1, pos2, frame.get_lattice(), self.bins, self.r_max, self.dr
            )

            if s1 == s2:
                n1 = len(pos1)
                norm = number_density * n1 / total_atoms
            else:
                n1 = len(pos1)
                n2 = len(pos2)
                norm = number_density * n2 / total_atoms

            shell_volumes = 4.0 * np.pi * self.r**2 * self.dr
            # Avoid division by zero for the first bin
            shell_volumes[0] = 4.0 / 3.0 * np.pi * (self.dr / 2.0) ** 3

            g_r = hist / (shell_volumes * norm * len(pos1))
            gr[pair] = g_r

        self.gr = gr

    @staticmethod
    @njit(parallel=True, nogil=True, cache=True, fastmath=True)
    def _calculate_histogram(
        pos1: np.ndarray,
        pos2: np.ndarray,
        lattice: np.ndarray,
        bins: int,
        r_max: float,
        dr: float,
    ) -> np.ndarray:
        hist = np.zeros(bins)
        for i in prange(len(pos1)):
            for j in range(len(pos2)):
                if np.array_equal(pos1[i], pos2[j]):
                    continue
                dist = calculate_pbc_distance(pos1[i], pos2[j], lattice)
                if dist < r_max:
                    bin_index = int(dist / dr)
                    hist[bin_index] += 1
        return hist
