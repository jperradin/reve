from reve.analysis.analyzer import Analyzer
from reve.io.reader import file_reader
from . import io
from . import core
from . import utils
from .settings.settings import Settings

from tqdm import tqdm

def main(settings: Settings):
    # Settings for files importation
    project_name = settings.project_name.value
    file_locations = settings.file_locations.value
    file_format = settings.file_format.value
    
    file_finder = io.FileFinder(file_locations, file_format)
    print(file_finder)
    files = file_finder.files
    
    file_reader = io.FileReader(files, file_format)
    print(file_reader)

    # Settings for analyzers
    analyzers = []
    for analyzer_type in settings.analyzers.value:
        analyzers.append(Analyzer.create(analyzer_type, settings))
    
    # Read first frame (TEST)
    lattices = file_reader.read_lattice(files[0], settings)
    atoms_frame = file_reader.read_atoms(files[0], settings)

    system = core.System(project_name, settings)
    system.set_current_frame(atoms_frame)
    print(system)

    # progress_bar = tqdm(
    #     range(number_of_frames), desc="Reading frames",
    #      unit="frame", leave=True, disable=quiet)
    
    # for _ in progress_bar:
    #     atoms_frame = file_reader.read_atoms(files[0], settings)

    #     # Analyze frame
    #     for analyzer in analyzers:
    #         analyzer.analyze(atoms_frame)

    #     # Update current frame at the end of the iteration
    #     settings.current_frame.value += 1
    