import abc
from typing import Dict, Type, List

class FileValidator(abc.ABC):
    """Abstract base class for file format validators."""
    
    @abc.abstractmethod
    def validate_file(self, file_path: str) -> None:
        """Validate if the file is in the correct format.
        
        Args:
            file_path: Path to the file to validate
            
        Raises:
            ValueError: If the file is not in the correct format
        """
        pass

    @abc.abstractmethod
    def validate_comment_line(self, line: str) -> bool:
        """ Validate if the comment line corresponds to the correct format."""
        pass


class XYZFileValidator(FileValidator):
    """Validator for extended XYZ file format."""
    
    def validate_file(self, file_path: str) -> None:
        """Validate if the file is in the extended XYZ format.
        
        Args:
            file_path: Path to the file to validate
            
        Raises:
            ValueError: If the file is not in the correct format
        """
        with open(file_path, "r") as f:
            # line 1 is number of atoms
            # line 2 is comment line and must contain the Lattice information:
            # ... Lattice="lx 0 0 0 ly 0 0 0 lz" ...
            headers = f.readlines()[:2]
            n = headers[0].strip()
            if not n.isdigit():
                raise ValueError(f"{file_path} : line 0 (number of atoms) is not an integer, exiting ... ")
            comment = headers[1].strip()
            if "Lattice" not in comment:
                raise ValueError(f"{file_path} : line 1 (comment) is not in the correct format, missing lattice informations, exiting ... ")
            
            # next lines are atomic positions ...
    def validate_comment_line(self, line: str) -> bool:
        """ Validate if the comment line corresponds to the extended XYZ format."""
        if "Lattice" not in line:
            return False
        return True

class LAMMPSFileValidator(FileValidator):
    """Validator for LAMMPS data file format."""
    
    def validate_file(self, file_path: str) -> None:
        """Validate if the file is in the LAMMPS data format.
        
        Args:
            file_path: Path to the file to validate
            
        Raises:
            ValueError: If the file is not in the correct format
        """
        with open(file_path, "r") as f:
            lines = f.readlines()
            
            # Check line 1: should contain number of atoms
            if len(lines) < 1:
                raise ValueError(f"{file_path}: File is empty, exiting...")
                
            # First line should contain number of atoms
            first_line = lines[0].strip()
            if not any(char.isdigit() for char in first_line):
                raise ValueError(f"{file_path}: Line 1 does not contain the number of atoms, exiting...")
                
            # Check line 2: should contain "Atoms. Timestep:"
            if len(lines) < 2:
                raise ValueError(f"{file_path}: File has only one line, missing header information, exiting...")
                
            second_line = lines[1].strip()
            if "Atoms. Timestep:" not in second_line:
                raise ValueError(f"{file_path}: Line 2 does not contain 'Atoms. Timestep:', exiting...")

    def validate_comment_line(self, line: str) -> bool:
        """ Validate if the comment line corresponds to the LAMMPS format."""
        if "Atoms. Timestep:" not in line:
            return False
        return True

class FileValidatorFactory:
    """Factory for creating file validators based on file format."""
    
    _validators: Dict[str, Type[FileValidator]] = {
        "xyz": XYZFileValidator,
        "lammps": LAMMPSFileValidator
    }
    
    @classmethod
    def get_validator(cls, file_format: str) -> FileValidator:
        """Get a validator for the specified file format.
        
        Args:
            file_format: The file format to get a validator for
            
        Returns:
            A validator for the specified file format
            
        Raises:
            ValueError: If the file format is not supported
        """
        validator_class = cls._validators.get(file_format)
        if validator_class is None:
            raise ValueError(f"File format {file_format} is not supported, exiting...")
        
        return validator_class()
    
    @classmethod
    def supported_formats(cls) -> List[str]:
        """Get a list of supported file formats.
        
        Returns:
            A list of supported file formats
        """
        return list(cls._validators.keys())
