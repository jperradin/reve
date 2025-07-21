import numpy as np
import os
from itertools import combinations
from tqdm import tqdm
from typing import Dict, Optional, List

from .base_analyzer import BaseAnalyzer
from ...core.frame import Frame
from ...core.node import Node
from ...config.settings import Settings
from ...utils.geometry import calculate_pbc_angle


class BondAngularDistributionAnalyzer(BaseAnalyzer):
    """
    Computes the bond angular distribution for specified triplets in the system.
    """

    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self.angles: Optional[Dict[str, np.ndarray]] = None
        self.angle_data: list = []
        self.frame_count: int = 0
        self.angle_max: float = 180.0
        self.bins: int = 800
        self.d_angle: float = self.angle_max / self.bins

        # NOTE: This is temporary, this will go to the settings
        # Add any other triplets you want to analyze here
        self.triplets_to_analyze = ["O-Si-O", "Si-O-Si", "O-O-O", "Si-Si-Si"]

    def analyze(self, frame: Frame) -> None:
        angles_for_frame = {triplet: [] for triplet in self.triplets_to_analyze}

        progress_bar_kwargs = {
            "disable": not self._settings.verbose,
            "leave": False,
            "ncols": os.get_terminal_size().columns,
            "colour": "blue",
        }

        progress_bar = tqdm(
            enumerate(frame.nodes),
            desc="Calculating angles ...",
            unit="atom",
            initial=0,
            total=len(frame.nodes),
            **progress_bar_kwargs,
        )
        for i, node in progress_bar:
            if len(node.neighbors) < 2:
                continue

            for neighbor1, neighbor2 in combinations(node.neighbors, 2):
                if neighbor1.node_id != neighbor2.node_id:
                    for triplet_type in self.triplets_to_analyze:
                        neighbor_atom1, central_atom, neighbor_atom2 = (
                            triplet_type.split("-")
                        )

                        # Check if the triplet matches the one to be analyzed
                        if node.symbol == central_atom and (
                            (
                                neighbor1.symbol == neighbor_atom1
                                and neighbor2.symbol == neighbor_atom2
                            )
                            or (
                                neighbor1.symbol == neighbor_atom2
                                and neighbor2.symbol == neighbor_atom1
                            )
                        ):
                            angle = calculate_pbc_angle(
                                neighbor1.position,
                                node.position,
                                neighbor2.position,
                                frame.lattice,
                            )
                            angles_for_frame[triplet_type].append(angle)

        # Histogram the angles for the current frame
        histograms = {}
        for triplet, angle_list in angles_for_frame.items():
            hist, _ = np.histogram(
                angle_list, bins=self.bins, range=(0, self.angle_max)
            )
            histograms[triplet] = hist

        if histograms:
            self.angle_data.append(histograms)
            self.frame_count += 1

    def finalize(self) -> None:
        if self.frame_count == 0:
            return

        # Sum the histograms from all frames
        summed_hist = {
            triplet: np.zeros(self.bins) for triplet in self.triplets_to_analyze
        }
        for frame_hist in self.angle_data:
            for triplet, hist in frame_hist.items():
                summed_hist[triplet] += hist

        # Normalize the distributions
        self.angles = {}
        for triplet, hist in summed_hist.items():
            total_angles = np.sum(hist)
            if total_angles > 0:
                self.angles[triplet] = hist / (
                    total_angles * self.d_angle
                )  # Normalize by integral
            else:
                self.angles[triplet] = np.zeros_like(hist)

    def get_result(self) -> Dict[str, np.ndarray]:
        return self.angles if self.angles is not None else {}

    def print_to_file(self) -> None:
        self.finalize()
        if self.angles is None:
            return
        output_path = os.path.join(
            self._settings.export_directory, "bond_angular_distribution.dat"
        )
        angle_bins = (
            np.linspace(0, self.angle_max, self.bins, endpoint=False)
            + self.d_angle / 2.0
        )

        with open(output_path, "w") as f:
            header = "# Angle\t" + "\t".join(self.angles.keys()) + "\n"
            f.write(header)
            for i in range(len(angle_bins)):
                line = f"{angle_bins[i]:.2f}"
                for triplet in self.angles.keys():
                    line += f"\t{self.angles[triplet][i]:.5f}"
                f.write(line + "\n")
