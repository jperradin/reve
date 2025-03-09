from typing import List, Dict, Type
from ..parser.file_validator import FileValidator, FileValidatorFactory
from .readers import Reader, XYZFileReader, LAMMPSFileReader
from ...settings.settings import Settings   
from ...core.frame import Frame

class FileReader:
    def __init__(self, files: List[str], file_format:str) -> None:
        self.files: List[str] = files
        self.file_format: str = file_format
        self.validator: FileValidator = FileValidatorFactory.get_validator(file_format)
        self.reader: Reader = FileReaderFactory.get_reader(file_format, files, self.validator)

    def read_atoms(self, file_path: str, settings: Settings) -> Frame:
        return self.reader.read_atoms(file_path, settings)

    def read_lattice(self, file_path: str, settings: Settings) -> Dict:
        return self.reader.read_lattice(file_path, settings)

    def __str__(self) -> str:
        return f"FileReader(files='{self.files}', file_format='{self.file_format}')"


class FileReaderFactory:
    """Factory for creating file readers based on file format."""

    _readers: Dict[str, Type[Reader]] = {
        "xyz": XYZFileReader,
        "lammps": LAMMPSFileReader
    }

    @classmethod
    def get_reader(cls, file_format: str, files: List[str], validator: FileValidator) -> Reader:
        """Get a reader for the specified file format.
        
        Args:
            file_format: The file format to get a reader for
            
        Returns:
            A reader for the specified file format
            
        Raises:
            ValueError: If the file format is not supported
        """
        reader_class = cls._readers.get(file_format)
        if reader_class is None:
            raise ValueError(f"File format {file_format} is not supported, exiting ... ")
        
        return reader_class(files, validator)