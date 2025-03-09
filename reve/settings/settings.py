from .parameter import Setting
from typing import Tuple

class Settings:
    def __init__(self) -> None:
        self.project_name: Setting = Setting("project_name", "default")
        self.export_directory: Setting = Setting("export_directory", "export")
        self._output_directory: str = ""
        self.file_locations: Setting = Setting("file_locations", "./")
        self.file_format: Setting = Setting("file_format", "xyz")
        self.number_of_atoms: Setting = Setting("number_of_atoms", 0)
        self.number_of_frames: Setting = Setting("number_of_frames", 0)
        self.range_of_frames: Setting = Setting("range_of_frames", (0, 0))
        self.current_frame: Setting = Setting("current_frame", 0)
        self.timestep: Setting = Setting("timestep", 0.0016)
        self.temperature: Setting = Setting("temperature", 300)
        self.pressure: Setting = Setting("pressure", 0)
        self.volume: Setting = Setting("volume", 0)
        self.quiet: Setting = Setting("quiet", False)
        self.overwrite_results: Setting = Setting("overwrite_results", False)
        self.analyzers: Setting = Setting("analyzers", [])

        # Composition settings
        self.composition: Setting = Setting("composition", "SiO2") # Default composition

        # Analyzers settings
        self.pdf_settings: Setting = Setting("pdf_settings", {"r_max": 10.0, "bin_width": 0.05})
        # self.bond_angular_settings: Setting = Setting("bond_angular_settings", {})
        # self.neutron_structure_factor_settings: Setting = Setting("neutron_structure_factor_settings", {})