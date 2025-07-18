from typing import List, Dict
from ...core.frame import Frame
from ...core.system import Settings
from .base_analyzer import BaseAnalyzer
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os
from datetime import datetime


class NeutronStructureFactorAnalyzer(BaseAnalyzer):
    """
    Computes the neutron structure factor for each pair of atoms
    """

    def __init__(self, settings: Settings) -> None:
        """Initializes the analyzer."""
        super().__init__(settings)

    def analyze(self, frame: Frame) -> None:
        return super().analyze(frame)

    def finalize(self) -> None:
        return super().finalize()

    def get_result(self) -> Dict[str, float]:
        return super().get_result()

    def print_to_file(self) -> None:
        return super().print_to_file()

    def _initialize_q_vectors(
        self, q_max: float = 10.0, q_step_scale: float = 6.0
    ) -> Dict:
        """
        Generate reciprocal space q-vectors.

        Args:
            q_max (float): Maximum magnitude for q-vectors. (less than half of the lattice dimension)
            q_step_scale (float): Scaling factor to determine the density of q-points.

        Returns:
            Dict: A dictionary containing, qx, qy, qz meshes, and their norms.
        """
        pass
