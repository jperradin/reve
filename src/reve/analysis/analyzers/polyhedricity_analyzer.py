import numpy as np
import os
from tqdm import tqdm
from typing import Dict, Optional, Tuple, List
from numba_progress import ProgressBar

from .base_analyzer import BaseAnalyzer
from ...core.frame import Frame
from ...core.node import Node
from ...config.settings import Settings
from ...utils.geometry import (
    calculate_octahedricity,
    calculate_pbc_dot_distances_combinations,
    calculate_square_based_pyramid,
    calculate_tetrahedricity,
    calculate_triangular_bipyramid,
)


class PolyhedricityAnalyzer(BaseAnalyzer):
    """Polyhedricity analyzer."""

    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self.polyhedricity: Optional[Dict[str, np.ndarray]] = None
        self.proportion: Optional[Dict[str, float]] = None
        self._atoms_data: Optional[Dict[str, np.ndarray]] = None
        self.central_nodes: List[Node] = []
        self.poly_data: List[Dict[str, np.ndarray]] = []
        self.counts = {
            "4_fold": 0,
            "5_fold": 0,
            "6_fold": 0,
            "5_fold_sbp": 0,
            "5_fold_tbp": 0,
        }

        # POLYAnalysisSettings
        self.central_species: str = (
            self._settings.analysis.poly_settings.central_species
            if self._settings.analysis.poly_settings is not None
            else "Si"
        )

        self.vertices_species: str = (
            self._settings.analysis.poly_settings.vertices_species
            if self._settings.analysis.poly_settings is not None
            else "O"
        )
        self.max_c: float = (
            self._settings.analysis.poly_settings.max_c
            if self._settings.analysis.poly_settings is not None
            else 0.5
        )
        self._bins: Optional[np.ndarray] = None
        self._dbin: Optional[float] = None
        self._mid: Optional[np.ndarray] = None

    def _initialize_arrays(self) -> None:
        """Pre-compute arrays that don't change between calculations."""
        if self._bins is None:
            self._bins = np.linspace(0, self.max_c, 1000)
            self._dbin = self.max_c / 1000
            self._mid = (self._bins[:-1] + self._bins[1:]) / 2
            # Tetrahedra
            self._hist_4_fold = np.zeros(len(self._bins))
            # all 5 folds considered as square based pyramid
            self._hist_SBP_pentahedricity = np.zeros(len(self._bins))
            # all 5 folds considered as triangular bipyramid
            self._hist_TBP_pentahedricity = np.zeros(len(self._bins))
            # minimum between square based pyramid vs triangular bipyramid
            self._hist_5_fold = np.zeros(len(self._bins))
            self._hist_sbp_pentahedricity = np.zeros(len(self._bins))
            self._hist_tbp_pentahedricity = np.zeros(len(self._bins))
            # Octahedra
            self._hist_6_fold = np.zeros(len(self._bins))

    def analyze(self, frame: Frame) -> None:
        self._initialize_arrays()
        self._atoms_data = frame.nodes_data.wrapped_positions
        self.central_nodes = frame.get_nodes_by_element(self.central_species)
        lattice = frame.get_lattice()
        N = len(self.central_nodes)
        self.polyhedricity = {}

        progress_bar_kwargs = {
            "disable": not self._settings.verbose,
            "leave": False,
            "ncols": os.get_terminal_size().columns,
            "colour": "blue",
        }

        progress_bar = tqdm(
            enumerate(self.central_nodes),
            desc="Calculating polyhedricity ...",
            unit="atom",
            initial=0,
            total=N,
            **progress_bar_kwargs,
        )
        for i, node in progress_bar:
            pos_batch = node.get_neighbors_positions_by_element(self.vertices_species)
            distances = calculate_pbc_dot_distances_combinations(pos_batch, lattice)
            distances.sort()

            polyhedricity = self._calculate_polyhedricity(node, distances)

            if polyhedricity is None:
                continue

            m, n = polyhedricity

            if int(m / self._dbin) + 1 > max(self._bins):
                continue

            if node.coordination == 4:
                self._hist_4_fold[int(m / self._dbin) + 1] += 1
                self.counts["4_fold"] += 1
            if node.coordination == 5:
                self.counts["5_fold"] += 1
                self._hist_SBP_pentahedricity[int(m / self._dbin) + 1] += 1
                self._hist_TBP_pentahedricity[int(n / self._dbin) + 1] += 1
                if m < n:
                    self.counts["5_fold_sbp"] += 1
                    self._hist_5_fold[int(m / self._dbin) + 1] += 1
                    self._hist_sbp_pentahedricity[int(m / self._dbin) + 1] += 1
                else:
                    self.counts["5_fold_tbp"] += 1
                    self._hist_5_fold[int(n / self._dbin) + 1] += 1
                    self._hist_tbp_pentahedricity[int(n / self._dbin) + 1] += 1
            if node.coordination == 6:
                self.counts["6_fold"] += 1
                self._hist_6_fold[int(m / self._dbin) + 1] += 1

        self.polyhedricity["4_fold"] = self._hist_4_fold
        self.polyhedricity["5_fold"] = self._hist_5_fold
        self.polyhedricity["6_fold"] = self._hist_6_fold
        self.polyhedricity["5_fold_sbp"] = self._hist_sbp_pentahedricity
        self.polyhedricity["5_fold_tbp"] = self._hist_tbp_pentahedricity
        self.polyhedricity["5_fold_SBP"] = self._hist_SBP_pentahedricity
        self.polyhedricity["5_fold_TBP"] = self._hist_TBP_pentahedricity

        self.poly_data.append(dict(self.polyhedricity))
        self.frame_processed_count += 1

    def finalize(self) -> None:
        self.counts["5_fold_SBP"] = self.counts["5_fold"]
        self.counts["5_fold_TBP"] = self.counts["5_fold"]

        if self.poly_data:
            keys = list(self.poly_data[0].keys())
            self.polyhedricity = {}
            self.proportion = {}

            self.polyhedricity["bins"] = self._bins

            for key in keys:
                values = np.array([frame_data[key] for frame_data in self.poly_data])
                # Avoid zero division
                if self.counts[key] == 0:
                    self.counts[key] = 1
                self.polyhedricity[key] = np.mean(values, axis=0) / self.counts[key]
                self.proportion[key] = self.counts[key] / (
                    len(self.central_nodes) * self.frame_processed_count
                )

    def get_result(self) -> Dict[str, float]:
        return self.proportion

    def print_to_file(self) -> None:
        self.finalize()
        if self.polyhedricity is None:
            return

        output_path1 = os.path.join(
            self._settings.export_directory, "polyhedricity_proportions.dat"
        )
        output_path2 = os.path.join(
            self._settings.export_directory, "polyhedricity_histograms.dat"
        )

        keys = list(self.proportion.keys())
        data = np.column_stack([self.proportion[k] for k in keys])

        np.savetxt(
            output_path1,
            data,
            header=f"{keys}",
            delimiter="\t",
            fmt="%.5f",
            comments="# ",
        )

        keys = list(self.polyhedricity.keys())
        data = np.column_stack([self.polyhedricity[k] for k in keys])

        np.savetxt(
            output_path2,
            data,
            header=f"{keys}",
            delimiter="\t",
            fmt="%.5f",
            comments="# ",
        )

    def _calculate_polyhedricity(
        self, node: Node, distances: np.ndarray
    ) -> Tuple[float, float]:
        # Possible errors
        if node.coordination == 4 and len(distances) != 6:
            raise ValueError(
                f"4-fold units should have 6 distances, got {len(distances)}"
            )
        if node.coordination == 5 and len(distances) != 10:
            raise ValueError(
                f"5-fold units should have 10 distances, got {len(distances)}"
            )
        if node.coordination == 6 and len(distances) != 15:
            raise ValueError(
                f"6-fold units should have 15 distances, got {len(distances)}"
            )

        if node.coordination == 4 and len(distances) == 6:
            # Calculate tetrahedricity
            tetrahedricity = calculate_tetrahedricity(distances)
            return tetrahedricity, -1.0
        if node.coordination == 5 and len(distances) == 10:
            # Calculate square based pyramid polyhedricity
            sbp_polyhedricity = calculate_square_based_pyramid(distances)
            # Calculate triangular bipyramid polyhedricity
            tbp_polyhedricity = calculate_triangular_bipyramid(distances)
            return sbp_polyhedricity, tbp_polyhedricity
        if node.coordination == 6 and len(distances) == 15:
            # Calculate triangular bipyramid polyhedricity
            octahedricity = calculate_octahedricity(distances)
            return octahedricity, -1.0
