import numpy as np
from typing import Dict, Any, List, Tuple
from .analyzer import Analyzer
from ..core.frame import Frame
from ..settings.settings import Settings

@Analyzer.register("pair_distribution_function")
class PairDistributionFunctionAnalyzer(Analyzer):
    """ Analyzer for calculating the pair distribution function g(r). """

    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        # Initialize default parameter for g(r) calculation 
        # if not specified in settings
        self.r_max = settings.pdf_settings.r_max.value if hasattr(settings.pdf_settings, "r_max") else 10.0
        self.bin_width = settings.pdf_settings.bin_width.value if hasattr(settings.pdf_settings, "bin_width") else 0.05
        self.n_bins = int(self.r_max / self.bin_width)
        self.histogram = np.zeros(self.n_bins)
        self.frame_count = 0

    def analyze(self, frame: Frame) -> Dict[str, Any]:
        """ Calculate the g(r) for a given frame. """
        positions = frame.positions
        lattice = frame.lattice

        # Use numba-optimized distance calculation
        # This is a placeholder for the actual implementation
        # which would use efficient distance calculation with PBC
        
        # Update histogram
        # Implementation would go here
        
        return {"frame_id": frame.id, "status": "processed"}

    def finalize(self) -> Dict[str, Any]:
        """ Normalize histogram to get final g(r) """
        if self.frame_count == 0:
            return {"error": "No frames processed"}
            
        # Normalize histogram
        # Implementation would go here
        
        # Calculate r values for each bin
        r_values = np.linspace(0, self.r_max, self.n_bins, endpoint=False) + self.bin_width / 2

        # Return results
        results = {
            "r_values": r_values,
            "g_r": self.histogram, # keep this for now
            "metadata": {
                "r_max": self.r_max,
                "bin_width": self.bin_width,
                "number_of_frames": self.frame_count
            }
        }

        return results