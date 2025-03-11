from tqdm import tqdm

from .settings import Settings
from .io.reader_factory import ReaderFactory
from .core.system import System
from .utils import generate_color_gradient

def main(settings: Settings):
    """
    Main function to test the package.
    """

    if not settings.quiet:
        print(settings)
    reader = ReaderFactory().get_reader(settings.file_location)
    reader.set_quiet(settings.quiet)

    system = System(reader, settings)

    if settings.range_of_frames[1] == -1:
        total = system.get_num_frames()
    else:
        total = settings.range_of_frames[1] - settings.range_of_frames[0]

    progress_bar = tqdm(system.iter_frames(), desc="Parsing frames", unit="frame", disable=settings.quiet, leave=False, total=total)
    for _ in progress_bar:
        pass
    
