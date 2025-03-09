import pathlib
import glob
from typing import List
from .file_validator import FileValidatorFactory, FileValidator

class FileFinder:
    """ 
    Finds files in a directory with a specific format.
    Supported formats: xyz, lammps

    #### Attributes:
    ---------------
        - path (pathlib.Path) : The path to the directory containing the trajectory files.
        - files (List[str]) : The list of trajectory files.
        - file_format (str) : The format of the trajectory files.
        - validator (FileValidator) : The validator for the file format.
    """

    def __init__(self, path: str, file_format: str) -> None:
        """ Initializes the FileFinder object. 
        
        #### Parameters:
        ---------------
            - path (str) : The path to the directory containing the trajectory files.
            - file_format (str) : The format of the trajectory files.
        """
        self.path: pathlib.Path = pathlib.Path(path)
        self.file_format: str = file_format
        self.validator: FileValidator = None
        self.files: List[str] = []
        self.number_of_files: int = 0
        
        if not self.path.is_dir():
            raise ValueError(f"Path {self.path} is not a directory, exiting ... ")

        # Check if the file format is supported
        if file_format not in FileValidatorFactory.supported_formats():
            raise ValueError(f"File format {file_format} is not supported, exiting ... ")
        
        # Get the appropriate validator for the file format
        self.validator = FileValidatorFactory.get_validator(file_format)
        
        # Find files with the given format
        # For LAMMPS files, we need to use a different extension
        if file_format == "lammps":
            self.files: List[str] = glob.glob(str(self.path / "*.data"))
        else:
            self.files: List[str] = glob.glob(str(self.path / f"*.{file_format}"))

        if len(self.files) == 0:
            raise ValueError("No trajectory files found in the directory, exiting ... ")

        # Validate files at initialization
        for file in self.files:
            self._validate_trajectory_file_format(file)

    def _validate_trajectory_file_format(self, file: str) -> None:
        """ Reads the header of the trajectory file to check if the format is correct"""
        self.validator.validate_file(file)
        self.number_of_files += 1

    def __str__(self) -> str:
        return f"FileFinder(path='{self.path}', file_format='{self.file_format}', number_of_files={self.number_of_files})"
