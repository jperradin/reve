from typing import Optional
from reve.io.base_reader import BaseReader
from reve.io.xyz_reader import XYZReader


class ReaderFactory:
    """Factory for creating file readers based on file type."""

    def __init__(self):
        self._readers = {}
        self.register_reader(XYZReader())
        # self.register_reader(LAMMPSReader())  # Register other readers here

    def register_reader(self, reader: BaseReader):
        """Registers a new reader instance."""
        # Use a dummy filename with the correct extension to determine support
        for ext in ['.xyz', '.lammpstrj', '.other']: #add your extensions here.
            if reader.supports_file_format(f'dummy{ext}'):
                self._readers[ext] = reader
                break


    def get_reader(self, filename: str) -> Optional[BaseReader]:
        """Returns the appropriate reader for a given file."""
        for extension, reader in self._readers.items():
            if filename.lower().endswith(extension):
                return reader
        return None