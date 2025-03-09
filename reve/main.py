from .settings import Settings
from .io.reader_factory import ReaderFactory


def main(settings: Settings):
    """
    Main function to test the package.
    """

    if not settings.quiet:
        print(settings)
    reader = ReaderFactory().get_reader(settings.file_location)
    frames, metadata = reader.read(settings.file_location)

