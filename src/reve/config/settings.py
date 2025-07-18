import os
from dataclasses import dataclass, field
import numpy as np
from typing import Tuple, Optional, List


@dataclass
class Cutoff:
    """
    Cutoff that contains all the cutoffs.

    Attributes:
    """

    type1: str
    type2: str
    distance: float

    def __str__(self) -> str:
        max_len = 5
        diff = max_len - len(self.type1) - 1 - len(self.type2)
        return f"{self.type1}-{self.type2}{' ' * diff} : distance = {self.distance}"

    def get_distance(self) -> float:
        return self.distance


@dataclass
class GeneralSettings:
    """
    General settings that contains all the general settings.

    Attributes:
    """

    project_name: str = "Project"  # Name of the project
    export_directory: str = "exports"  # Directory to export results
    file_location: str = ""  # Path to the trajectory file
    range_of_frames: Tuple[int, int] = (
        0,
        -1,
    )  # Range of frames to process (0 to -1 = all frames)
    apply_pbc: bool = False  # Whether to apply periodic boundary conditions
    verbose: bool = (
        False  # Whether to print settings, progress bars and other information
    )
    save_logs: bool = False  # Whether to save logs
    save_performance: bool = False  # Whether to save performance
    cutoffs: List[Cutoff] = field(default_factory=lambda: [])  # Cutoffs for distance
    coordination_mode: str = "all_types"


@dataclass
class AnalysisSettings:
    """
    Analysis settings that contains all the analyzer settings.

    Attributes:
    """

    overwrite: bool = True  # Whether to overwrite the existing file, if False, appends results to the file
    with_all: bool = False  # Whether to calculate all the properties
    with_neutron_structure_factor: bool = (
        False  # Whether to calculate the neutron structure factor
    )

    def get_analyzers(self) -> List[str]:
        analyzers = []
        if self.with_neutron_structure_factor:
            analyzers.append("NeutronStructureFactorAnalyzer")
        if self.with_all:
            analyzers.append("NeutronStructureFactorAnalyzer")
        return analyzers

    def __str__(self) -> str:
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                if not self.with_all and key == "with_all":
                    continue
                elif (
                    not self.with_neutron_structure_factor
                    and key == "with_neutron_structure_factor"
                ):
                    continue
                lines.append(f"\t\t|- {key}: {value}")
        output = """
        Analysis Settings:
        -----------------
{}
        """.format("\n".join(lines))
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
        apply_pbc (bool): Whether to apply periodic boundary conditions.
    """

    apply_custom_lattice: bool = False
    custom_lattice: np.ndarray = field(
        default_factory=lambda: np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        )
    )
    get_lattice_from_file: bool = False
    lattice_file_location: str = "./"
    apply_lattice_to_all_frames: bool = True

    def __str__(self) -> str:
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                if not self.apply_custom_lattice and key == "apply_custom_lattice":
                    lines.append(f"\t\t|- {key}: {value}")
                    break
                elif key == "custom_lattice":
                    line1 = f"\t\t|- {key}:"
                    lx = np.array2string(
                        value[0],
                        separator=", ",
                        formatter={"float_kind": lambda x: f"{x}"},
                    )
                    ly = np.array2string(
                        value[1],
                        separator=", ",
                        formatter={"float_kind": lambda x: f"{x}"},
                    )
                    lz = np.array2string(
                        value[2],
                        separator=", ",
                        formatter={"float_kind": lambda x: f"{x}"},
                    )
                    lines.append(
                        f"{line1}\n\t\t\tlx = {lx}\n\t\t\tly = {ly}\n\t\t\tlz = {lz}"
                    )
                else:
                    lines.append(f"\t\t|- {key}: {value}")
        output = """

        Lattice Settings:
        -----------------
{}
        """.format("\n".join(lines))
        return output


@dataclass
class Settings:
    """Settings for the Reve package and it is constructed using the SettingsBuilder."""

    project_name: str = "default"
    export_directory: str = "export"
    file_location: str = "./"
    range_of_frames: Tuple[int, int] = (0, -1)
    apply_pbc: bool = True
    verbose: bool = False
    save_logs: bool = False
    save_performance: bool = False
    coordination_mode: str = "all_types"
    general: GeneralSettings = field(default_factory=GeneralSettings)
    cutoffs: List[Cutoff] = field(default_factory=lambda: [])
    lattice: LatticeSettings = field(default_factory=LatticeSettings)
    analysis: AnalysisSettings = field(default_factory=AnalysisSettings)

    @property
    def output_directory(self) -> str:
        return os.path.join(self.export_directory, self.project_name)

    def get_cutoff(self, type1: str, type2: str) -> float | None:
        for cutoff in self.cutoffs:
            if cutoff.type1 == type1 and cutoff.type2 == type2:
                return cutoff.distance
            elif cutoff.type1 == type2 and cutoff.type2 == type1:
                return cutoff.distance
        return None

    def get_max_cutoff(self) -> float:
        return max(c.distance for c in self.cutoffs)

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
                if key == "general":
                    continue
                elif key == "lattice":
                    lines.append(f"\t{str(self.lattice)}")
                elif key == "analysis":
                    lines.append(f"\t{str(self.analysis)}")
                elif key == "cutoffs":
                    line1 = f"\t|- {key:}:"
                    for cutoff in value:
                        line1 += f"\n\t\t{str(cutoff)}"
                    lines.append(line1)
                else:
                    lines.append(f"\t|- {key}: {value}")

        output = """
        General Settings:
        ----------------
{}
        """.format("\n".join(lines))
        return output


class SettingsBuilder:
    def __init__(self):
        self._settings = Settings()  # Start with default settings

    def with_lattice(self, lattice: LatticeSettings):
        if not isinstance(lattice, LatticeSettings):
            raise ValueError(f"Invalid lattice settings: {lattice}")
        self._settings.lattice = lattice
        return self

    def with_general(self, general: GeneralSettings):
        if not isinstance(general, GeneralSettings):
            raise ValueError(f"Invalid general settings: {general}")
        if not general.project_name:
            raise ValueError(f"Invalid project name: {general.project_name}")
        if not general.export_directory:
            raise ValueError(f"Invalid export directory: {general.export_directory}")
        if not general.file_location:
            raise ValueError(f"Invalid file location: {general.file_location}")
        if not general.range_of_frames:
            raise ValueError(f"Invalid range of frames: {general.range_of_frames}")
        if not general.cutoffs:
            raise ValueError(f"Invalid cutoffs list: {general.cutoffs}")
        if general.apply_pbc is None:
            raise ValueError(f"Invalid apply pbc: {general.apply_pbc}")

        self._settings.project_name = general.project_name
        self._settings.export_directory = general.export_directory
        self._settings.file_location = general.file_location
        self._settings.range_of_frames = general.range_of_frames
        self._settings.apply_pbc = general.apply_pbc
        self._settings.cutoffs = general.cutoffs
        self._settings.coordination_mode = general.coordination_mode
        if general.verbose is not None:
            self._settings.verbose = general.verbose
        if general.save_logs is not None:
            self._settings.save_logs = general.save_logs
        if general.save_performance is not None:
            self._settings.save_performance = general.save_performance
        return self

    def with_analysis(self, analysis: AnalysisSettings):
        if not isinstance(analysis, AnalysisSettings):
            raise ValueError(f"Invalid analysis settings: {analysis}")
        self._settings.analysis = analysis
        return self

    def build(self) -> Settings:
        return self._settings


__all__ = [
    Settings,
    SettingsBuilder,
    AnalysisSettings,
    LatticeSettings,
    Cutoff,
]
