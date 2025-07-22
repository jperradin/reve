import os
from dataclasses import dataclass, field
import numpy as np
from typing import Tuple, Optional, List


@dataclass
class PDFAnalysisSettings:
    """Settings specific to the PairDistributionFunctionAnalyzer"""

    r_max: float = 10.0
    bins: int = 800
    pairs_to_calculate: List[str] = field(default_factory=lambda: [""])

    def __str__(self) -> str:
        line = "\t\t|- pdf_settings:\n"
        line += f"\t\t  |- r_max = {self.r_max}\n"
        line += f"\t\t  |- bins = {self.bins}\n"
        line += f"\t\t  |- pairs_to_calculate = \n\t\t    {self.pairs_to_calculate}"
        return line


@dataclass
class BADAnalysisSettings:
    """Settings specific to the BondAngularDistributionAnalyzer"""

    bins: int = 800
    triplets_to_calculate: List[str] = field(default_factory=lambda: [""])

    def __str__(self) -> str:
        line = "\t\t|- bad_settings:\n"
        line += f"\t\t  |- bins = {self.bins}\n"
        line += (
            f"\t\t  |- triplets_to_calculate = \n\t\t    {self.triplets_to_calculate}"
        )
        return line


@dataclass
class SQFFTAnalysisSettings:
    """Settings specific to the NeutronStructureFactorFFTAnalyzer"""

    grid_size_max: int = 512


@dataclass
class STRUNITSAnalysisSettings:
    """Settings specific to the StructuralUnitsAnalyzer"""

    units_to_calculate: List[str] = field(default_factory=lambda: [""])

    def __str__(self) -> str:
        line = "\t\t|- strunits_settings:\n"
        line += f"\t\t  |- units_to_calculate = \n\t\t    {self.units_to_calculate}"
        return line


@dataclass
class CONNAnalysisSettings:
    """Settings specific to the ConnectivityAnalyzer"""

    networking_species: str = ""
    bridging_species: str = ""

    def __str__(self) -> str:
        line = "\t\t|- connect_settings:\n"
        line += f"\t\t  |- networking_species = \n\t\t    {self.networking_species}"
        line += f"\t\t  |- bridging_species = \n\t\t    {self.bridging_species}"
        return line


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

    # Name of the project
    project_name: str = "Project"
    # Directory to export results
    export_directory: str = "exports"
    # Path to the trajectory file
    file_location: str = ""
    # Range of frames to process (0 to -1 = all frames)
    range_of_frames: Tuple[int, int] = (0, -1)
    # Whether to apply periodic boundary conditions
    apply_pbc: bool = False
    # Whether to wrap positions systematically
    wrap_position: bool = True
    # Whether to print settings, progress bars and other information
    verbose: bool = False
    # Whether to save logs
    save_logs: bool = False
    # Whether to save performance
    save_performance: bool = False
    # Cutoffs for distance
    cutoffs: List[Cutoff] = field(default_factory=lambda: [])
    # Coordination mode
    coordination_mode: str = "all_types"


@dataclass
class AnalysisSettings:
    """
    Analysis settings that contains all the analyzer settings.

    Attributes:
    """

    # Whether to overwrite the existing file, if False, appends results to the file
    overwrite: bool = True
    # Whether to calculate all the properties
    with_all: bool = False
    # Whether to calculate the pair distribution functions.
    with_pair_distribution_function: bool = False
    pdf_settings: Optional[PDFAnalysisSettings] = None
    # Whether to calculate the bond angular distribution
    with_bond_angular_distribution: bool = False
    bad_settings: Optional[BADAnalysisSettings] = None
    # Whether to calculate the neutron structure factor
    with_neutron_structure_factor: bool = False
    # Whether to calculate the neutron structure factor via fft
    with_neutron_structure_factor_fft: bool = False
    sqfft_settings: Optional[SQFFTAnalysisSettings] = None
    # Whether to calculate the structural units
    with_structural_units: bool = False
    strunits_settings: Optional[STRUNITSAnalysisSettings] = None
    # Whether to calculate the conncetivities
    with_connectivity: bool = False
    connect_settings: Optional[CONNAnalysisSettings] = None

    def get_analyzers(self) -> List[str]:
        analyzers = []
        if self.with_pair_distribution_function:
            analyzers.append("PairDistributionFunctionAnalyzer")
        if self.with_bond_angular_distribution:
            analyzers.append("BondAngularDistributionAnalyzer")
        if self.with_neutron_structure_factor:
            raise ValueError(
                "The NeutronStructureFactorAnalyzer is disable, use NeutronStructureFactorFFTAnalyzer instead"
            )
            # analyzers.append("NeutronStructureFactorAnalyzer")
        if self.with_neutron_structure_factor_fft:
            analyzers.append("NeutronStructureFactorFFTAnalyzer")
        if self.with_structural_units:
            analyzers.append("StructuralUnitsAnalyzer")
        if self.with_connectivity:
            analyzers.append("ConnectivityAnalyzer")
        if self.with_all:
            analyzers.append("PairDistributionFunctionAnalyzer")
            analyzers.append("BondAngularDistributionAnalyzer")
            analyzers.append("NeutronStructureFactorFFTAnalyzer")
            analyzers.append("StructuralUnitsAnalyzer")
            analyzers.append("ConnectivityAnalyzer")
            # analyzers.append("NeutronStructureFactorAnalyzer")
        return analyzers

    def __str__(self) -> str:
        lines = []
        for key, value in self.__dict__.items():
            hold = key, value
            if value is not None:
                if not self.with_all and key == "with_all":
                    continue
                if not self.with_pair_distribution_function and (
                    key == "with_pair_distribution_function" or key == "pdf_settings"
                ):
                    continue
                if not self.with_bond_angular_distribution and (
                    key == "with_bond_angular_distribution" or key == "bad_settings"
                ):
                    continue
                if (
                    not self.with_neutron_structure_factor
                    and key == "with_neutron_structure_factor"
                ):
                    continue
                if not self.with_neutron_structure_factor_fft and (
                    key == "with_neutron_structure_factor_fft"
                    or key == "sqfft_settings"
                ):
                    continue
                if not self.with_structural_units and (
                    key == "with_structural_units" or key == "strunits_settings"
                ):
                    continue
                if not self.with_connectivity and (
                    key == "with_connectivity" or key == "connect_settings"
                ):
                    continue
                if (
                    self.with_pair_distribution_function or self.with_all
                ) and key == "pdf_settings":
                    lines.append(str(self.pdf_settings))
                    continue
                if (
                    self.with_bond_angular_distribution or self.with_all
                ) and key == "bad_settings":
                    lines.append(str(self.bad_settings))
                    continue
                if (
                    self.with_structural_units or self.with_all
                ) and key == "strunits_settings":
                    lines.append(str(self.strunits_settings))
                    continue
                if (
                    self.with_connectivity or self.with_all
                ) and key == "connect_settings":
                    lines.append(str(self.connect_settings))
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
    wrap_position: bool = True
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
        if (
            analysis.with_bond_angular_distribution or analysis.with_all
        ) and analysis.bad_settings is None:
            # Create a default BADAnalysisSettings object if not created
            bad_settings = BADAnalysisSettings(
                triplets_to_calculate=["O-O-O", "O-Si-O", "Si-O-Si", "Si-Si-Si"]
            )
            analysis.bad_settings = bad_settings
        if (
            analysis.with_pair_distribution_function or analysis.with_all
        ) and analysis.pdf_settings is None:
            # Create a default PDFAnalysisSettings object if not created
            pdf_settings = PDFAnalysisSettings(
                pairs_to_calculate=["O-O", "O-Si", "Si-Si", "total"]
            )
            analysis.pdf_settings = pdf_settings
        if (
            analysis.with_neutron_structure_factor_fft
            and analysis.sqfft_settings is None
        ):
            # Create a SQFFTAnalysisSettings object if not create
            sqfft_settings = SQFFTAnalysisSettings()
            analysis.sqfft_settings = sqfft_settings
        if analysis.with_structural_units and analysis.strunits_settings is None:
            # Create a SQFFTAnalysisSettings object if not create
            strunits_settings = STRUNITSAnalysisSettings(
                units_to_calculate=[
                    "SiO4",
                    "SiO5",
                    "SiO6",
                    "SiO7",
                    "OSi1",
                    "OSi2",
                    "OSi3",
                    "OSi4",
                ]
            )
            analysis.strunits_settings = strunits_settings
        if analysis.with_connectivity and analysis.connect_settings is None:
            # Create a SQFFTAnalysisSettings object if not create
            connect_settings = CONNAnalysisSettings(
                networking_species="Si",
                bridging_species="O",
            )
            analysis.connect_settings = connect_settings

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
    GeneralSettings,
    PDFAnalysisSettings,
    BADAnalysisSettings,
    SQFFTAnalysisSettings,
    STRUNITSAnalysisSettings,
]
