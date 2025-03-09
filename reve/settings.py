import os
from dataclasses import dataclass, field
from typing import Tuple, List, Dict, Optional, override

@dataclass
class PDFSettings:
    r_max: float = 10.0
    bin_width: float = 0.05

@dataclass
class BDFSettings:  # Example
    angle_max: float = 180.0
    bin_width: float = 1.0

@dataclass
class AnalysisSettings:
    pdf: PDFSettings = field(default_factory=PDFSettings)
    bdf: BDFSettings = field(default_factory=BDFSettings)
    # Add other analysis settings here

@dataclass
class Settings:
    project_name: str = "default"
    export_directory: str = "export"
    file_location: str = "./"
    number_of_atoms: int = 0
    range_of_frames: Tuple[int, int] = (0, -1)
    timestep: Optional[float] = None
    temperature: Optional[float] = None
    pressure: Optional[float] = None
    volume: Optional[float] = None
    quiet: bool = False
    overwrite_results: bool = False
    analysis: AnalysisSettings = field(default_factory=AnalysisSettings)

    @property
    def output_directory(self) -> str:
        return os.path.join(self.export_directory, self.project_name)

    def set_range_of_frames(self, start: int, end: Optional[int] = None):
        if end is None:
            end = -1
        if start < 0:
            raise ValueError("Start frame cannot be negative")
        if end != -1 and start > end:
            raise ValueError("Start frame cannot be greater than end frame")
        self.range_of_frames = (start, end)

    def __str__(self) -> str:
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                if key == 'analysis':
                    for analysis_type, analysis_settings in value.__dict__.items():
                        for setting_key, setting_value in analysis_settings.__dict__.items():
                            if setting_value is not None:
                                lines.append(f"\t|- {analysis_type} {setting_key}: {setting_value}")
                else:
                    lines.append(f"\t|- {key}: {value}")
        output = '''
        Settings:
        ---------
{}
        '''.format('\n'.join(lines))
        return output

class SettingsBuilder:
    def __init__(self):
        self._settings = Settings()  # Start with default settings

    def with_project_name(self, name: str):
        self._settings.project_name = name
        return self

    def with_export_directory(self, directory: str):
        self._settings.export_directory = directory
        return self

    def with_file_location(self, location: str):
        self._settings.file_location = location
        return self

    def with_number_of_atoms(self, num_atoms: int):
        self._settings.number_of_atoms = num_atoms
        return self

    def with_range_of_frames(self, start: int, end: Optional[int] = None):
        self._settings.set_range_of_frames(start, end)
        return self

    def with_timestep(self, timestep: float):
        self._settings.timestep = timestep
        return self
    
    def with_temperature(self, temperature: float):
        self._settings.temperature = temperature
        return self

    def with_pressure(self, pressure: float):
        self._settings.pressure = pressure
        return self

    def with_volume(self, volume: float):
        self._settings.volume = volume
        return self

    def with_quiet(self, quiet: bool):
        self._settings.quiet = quiet
        return self

    def with_overwrite_results(self, overwrite: bool):
        self._settings.overwrite_results = overwrite
        return self

    def with_pdf_settings(self, r_max: Optional[float] = None, bin_width: Optional[float] = None):
        if r_max is not None:
            self._settings.analysis.pdf.r_max = r_max
        if bin_width is not None:
            self._settings.analysis.pdf.bin_width = bin_width
        return self

    def with_bdf_settings(self, angle_max: Optional[float] = None, bin_width: Optional[float] = None):
        if angle_max is not None:
            self._settings.analysis.bdf.angle_max = angle_max
        if bin_width is not None:
            self._settings.analysis.bdf.bin_width = bin_width
        return self

    # Add more 'with_' methods for other analysis settings as needed

    def build(self) -> Settings:
        return self._settings

