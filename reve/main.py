from .settings import Settings
from .io.reader_factory import ReaderFactory
from .core.system import System


def main(settings: Settings):
    """
    Main function to test the package.
    """

    if not settings.quiet:
        print(settings)
    reader = ReaderFactory().get_reader(settings.file_location)
    reader.set_quiet(settings.quiet)

    system = System(reader, settings)
    
    for frame in system.iter_frames():
        print(f"Frame {frame.frame_id} loaded.")
    
