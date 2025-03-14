
# Import necessary modules
from reve import SettingsBuilder, main
import reve.config.settings as c
import numpy as np

# Lattice settings
config_lattice = c.LatticeSettings(
    apply_custom_lattice=True,                      # If False, read lattice from trajectory file
    custom_lattice=np.array([[66.2574, 0.0, 0.0],     # Put a custom lattice here or None
                             [0.0, 66.2574, 0.0],     # lattice values MUST be float.
                             [0.0, 0.0, 66.2574]]), 
    apply_lattice_to_all_frames=True,               # Apply the same lattice to all frames
    get_lattice_from_file=False,                    # Get lattice from file (will read from lattice.dat if True)
    lattice_file_location="./",                     # Location of the lattice file
)

# PDF settings
config_pdf = c.PDFSettings(
    run=True,              # Set True to run PDF analysis
    r_max=12.0,            # Maximum radius
    bin_width=0.025        # Width of each bin
)

# BAD settings
config_bad = c.BADSettings(
    run=True,               # Set True to run BAD analysis
    angle_max=180.0,        # Maximum angle
    bin_width=1.0           # Width of each bin
)

# Analysis settings
config_analysis = c.AnalysisSettings(
    pdf=config_pdf,         # PDF settings
    bdf=config_bad          # BAD settings
)

# Path to the trajectory file
path = 'tests/input-10_frames.xyz'

# Settings builder
settings = (SettingsBuilder() \
    .with_project_name('test')         # Name of the project \
    .with_export_directory('export')   # Directory to export results \
    .with_file_location(path)          # Path to the trajectory file \
    .with_number_of_atoms(27216)       # Number of atoms in the trajectory \
    .with_range_of_frames(0, -1)       # Range of frames to process (0 to -1 = all frames) \
    .with_verbose(True)               # Whether to print settings, progress bars and other information (True = print) \
    .with_lattice(config_lattice)      # Lattice settings \
    .with_analysis(config_analysis)    # Analysis settings \

    .build()                           # Don't forget to build the settings object
)

# Run the main function to process the trajectory
main(settings) 
