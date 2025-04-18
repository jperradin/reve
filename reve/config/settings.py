import os
from dataclasses import dataclass, field
import numpy as np
from typing import Tuple, Optional

@dataclass
class PDFSettings:
    """ 
    PDF settings. 
    
    Attributes:
        r_max (float): Maximum radius. (should be less than half the box size to avoid overlap by periodicity)
        bin_width (float): Width of each bin.
    """
    run: bool = True
    r_max: float = 10.0
    bin_width: float = 0.05

@dataclass
class BADSettings:
    """ 
    BDF settings. 
    
    Attributes:
        angle_max (float): Maximum angle in degrees.
        bin_width (float): Width of each bin in degrees.
    """
    run: bool = True
    angle_max: float = 180.0
    bin_width: float = 1.0


@dataclass
class AnalysisSettings:
    """ 
    Analysis settings that contains all the analyzer settings. 
    
    Attributes:
        pdf (PDFSettings): PDF settings.
        bdf (BADSettings): BDF settings.
    """
    pdf: PDFSettings = field(default_factory=PDFSettings)
    bdf: BADSettings = field(default_factory=BADSettings)
    # TODO: Add other analysis settings here

    def __str__(self) -> str:
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                lines.append(f"\t\t|- {value}")
        output = '''
        Analysis Settings:
        -----------------
{}
        '''.format('\n'.join(lines))
        return output

@dataclass
class LatticeSettings:
    """ 
    Lattice settings. 
    
    TODO implement lattice fetcher from file
         implement the handling of lattice settings in the system

    Attributes:
        lattice_in_trajectory_file (bool): Whether the lattice is present in the trajectory file.
        lattice (np.ndarray): The lattice matrix.
        get_lattice_from_file (bool): Whether to get the lattice from a file.
        lattice_file_location (str): Location of the lattice file.
        apply_lattice_to_all_frames (bool): Whether to apply the lattice to all frames.
    """
    apply_custom_lattice: bool = False
    custom_lattice: np.ndarray = field(default_factory=lambda: np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]))
    get_lattice_from_file: bool = False
    lattice_file_location: str = "./"
    apply_lattice_to_all_frames: bool = True

    def __str__(self) -> str:
        if not self.apply_custom_lattice:
            return ""
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                if key == "custom_lattice":
                    line1 = f"\t\t|- {key}:"
                    lx = np.array2string(value[0], separator=', ', formatter={'float_kind': lambda x: f'{x}'})
                    ly = np.array2string(value[1], separator=', ', formatter={'float_kind': lambda x: f'{x}'})
                    lz = np.array2string(value[2], separator=', ', formatter={'float_kind': lambda x: f'{x}'})
                    lines.append(f"{line1}\n\t\t\tlx = {lx}\n\t\t\tly = {ly}\n\t\t\tlz = {lz}")
                else:
                    lines.append(f"\t\t|- {key}: {value}")
        output = '''

        Lattice Settings:
        -----------------
{}
        '''.format('\n'.join(lines))
        return output

@dataclass
class Settings:
    """ Settings for the Reve package and it is constructed using the SettingsBuilder. """
    project_name: str = "default"
    export_directory: str = "export"
    file_location: str = "./"
    number_of_atoms: int = 0
    range_of_frames: Tuple[int, int] = (0, -1)
    timestep: float = 0.0
    verbose: bool = False
    lattice: LatticeSettings = field(default_factory=LatticeSettings) # TODO: implement lattice fetcher from file
    # analysis: AnalysisSettings = field(default_factory=AnalysisSettings) # TODO: implement analyzer first

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
                if key == 'lattice' and not self.lattice.apply_custom_lattice:
                    continue
                elif key == 'lattice' and self.lattice.apply_custom_lattice:
                    lines.append(f"\t{str(self.lattice)}")
                elif key == 'analysis':
                    lines.append(f"\t{str(self.analysis)}")
                else:
                    lines.append(f"\t|- {key}: {value}")
        output = '''
        Global Settings:
        ----------------
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

    def with_lattice(self, lattice: LatticeSettings):
        self._settings.lattice = lattice
        return self

    def with_verbose(self, verbose: bool):
        self._settings.verbose = verbose
        return self

    def with_analysis(self, analysis: AnalysisSettings):
        self._settings.analysis = analysis
        return self

    def with_pdf(self, pdf: PDFSettings):
        self._settings.analysis.pdf = pdf
        return self

    def with_bdf(self, bdf: BADSettings):
        self._settings.analysis.bdf = bdf
        return self

    def build(self) -> Settings:
        return self._settings

__all__ = [
    Settings,
    SettingsBuilder,
    PDFSettings,
    BADSettings,
    AnalysisSettings
]