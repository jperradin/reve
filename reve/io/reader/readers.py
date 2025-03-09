from itertools import islice
from typing import Dict, List
from ..parser.file_validator import FileValidator
from ...settings.settings import Settings
import abc
import numpy as np
from ...core.atom import Atom
from ...core.frame import Frame

class Reader(abc.ABC):
    @abc.abstractmethod
    def read_lattice(self, file_path: str, settings: Settings) -> Dict:
        pass

    @abc.abstractmethod
    def read_atoms(self, file_path: str, settings: Settings) -> Dict:
        pass


class XYZFileReader(Reader):
    def __init__(self, files: List[str], validator: FileValidator) -> None:
        self.files: List[str] = files
        self.validator: FileValidator = validator
        self.comment_lines: List[str] = []
        self.lattices: Dict[int, List[np.float64]] = {}
        
    def parse_frames_header(self, file_path: str) -> None:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        f.close()
        self.comment_lines = [line for line in lines if 'Lattice' in line]
        del lines

    def read_lattice(self, file_path: str, settings: Settings) -> None:
        """ Read the lattice properties from the trajectory file.
        
        Parameters
        ----------
            - file_path (str) : Path to the trajectory file.
            - settings (Settings) : Settings for reading the file.
            
        Returns:
        --------
            - Dict : Dictionary containing the lattice properties.
        """
        # First make sure we have parsed the headers
        if not self.comment_lines:
            self.parse_frames_header(file_path)

        lattices = [line for line in self.comment_lines if "Lattice" in line]

        lattices_dict = {}

        for n, line in enumerate(lattices):
            current_lattice = line.split('\"')[1]
            lxx = float(current_lattice.split()[0])
            lxy = float(current_lattice.split()[1])
            lxz = float(current_lattice.split()[2])
            lyx = float(current_lattice.split()[3])
            lyy = float(current_lattice.split()[4])
            lyz = float(current_lattice.split()[5])
            lzx = float(current_lattice.split()[6])
            lzy = float(current_lattice.split()[7])
            lzz = float(current_lattice.split()[8])

            lattices_dict[n] = np.array([
                [lxx, lxy, lxz],
                [lyx, lyy, lyz],
                [lzx, lzy, lzz]
            ])

        self.lattices = lattices_dict
        return 

    def read_atoms(self, file_path: str, settings: Settings) -> Frame:
        """
        Read atomic positions from the trajectory file.
        
        Parameters
        ----------
            - file_path (str) : Path to the trajectory file.
            - settings (Settings) : Settings for reading the file.
            
        Returns:
        --------
            - Dict : Dictionary containing the atomic positions.
        """
        with open(file_path, 'r') as f:
            lines = f.readlines()
        f.close()

        # Select frame to read
        first_frame = settings.range_of_frames.value[0]
        last_frame = settings.range_of_frames.value[1]
        cf = settings.current_frame.value
        if cf < first_frame or cf > last_frame:
            raise ValueError(f"Frame {cf} is not in the range {first_frame}-{last_frame}.")
        
        n_atoms = settings.number_of_atoms.value
        n_header = 2
        start_line = cf * (n_header + n_atoms) # frames are 0-indexed

        # Save memory
        frame_lines = lines[start_line:start_line + (n_header + n_atoms)]
        header_lines = frame_lines[:n_header]
        if n_atoms != int(header_lines[0]):
            raise ValueError(f"Frame {cf} has {n_atoms} atoms, expected {int(header_lines[0])}.")

        atom_lines = frame_lines[n_header:]
        del lines
        
        positions = np.zeros((n_atoms, 3), dtype=np.float64)
        elements = np.zeros(n_atoms, dtype=np.string_)
        atoms_ids = np.zeros(n_atoms, dtype=np.int64)
        
        for n, line in enumerate(atom_lines):
            # line format: atom_type x y z
            parts = line.split()

            atom_type = parts[0]
            x = np.float64(parts[1])
            y = np.float64(parts[2])
            z = np.float64(parts[3])

            position = np.array([x, y, z])

            positions[n] = position
            elements[n] = atom_type
            atoms_ids[n] = n

        del frame_lines
        del header_lines
        del atom_lines

        return Frame(cf, self.lattices[cf], positions, elements, atoms_ids)

    def read_atoms(self, file_path: str, settings: Settings) -> Frame:
        """
        Read atomic positions from the trajectory file.
        
        Parameters
        ----------
            - file_path (str) : Path to the trajectory file.
            - settings (Settings) : Settings for reading the file.
            
        Returns:
        --------
            - Dict : Dictionary containing the atomic positions.
        """
        # Get frame info
        first_frame = settings.range_of_frames.value[0]
        last_frame = settings.range_of_frames.value[1]
        cf = settings.current_frame.value
        if cf < first_frame or cf > last_frame:
            raise ValueError(f"Frame {cf} is not in the range {first_frame}-{last_frame}.")
        
        n_atoms = settings.number_of_atoms.value
        N_HEADER = 2
        frame_size = n_atoms + N_HEADER
        start_line = cf * frame_size
        
        # Skip to the right frame and read only necessary lines
        with open(file_path, 'r') as f:
            # Skip to the start of the frame
            for _ in range(start_line):
                next(f)


            # Read the header to verify atom count
            n_atoms_in_file = int(next(f).strip())
            if n_atoms != n_atoms_in_file:
                raise ValueError(f"Frame {cf} has {n_atoms} atoms, expected {n_atoms_in_file}.")

            # Skip the comment line
            next(f)

            # Use numpy's loadtxt to efficiently read the data
            # Create a converter to handle the element column
            elements = np.zeros(n_atoms, dtype=np.string_)
            positions = np.zeros((n_atoms, 3), dtype=np.float64)

            for i, line in enumerate(islice(f, n_atoms)):
                parts = line.split()
                elements[i] = parts[0]
                positions[i] = np.array([np.float64(parts[1]), np.float64(parts[2]), np.float64(parts[3])])

        f.close()

        # Create atom ids
        atoms_ids = np.arange(n_atoms, dtype=np.int64)

        return Frame(cf, self.lattices[cf], positions, elements, atoms_ids)

        


class LAMMPSFileReader(Reader):
    def __init__(self, files: List[str], validator: FileValidator) -> None:
        self.files: List[str] = files
        self.validator: FileValidator = validator
        self.comment_lines: List[str] = []
        
    def read_lattice(self, file_path: str, settings: Settings) -> Dict:
        """
        Read the lattice properties from the trajectory file.
        
        Parameters
        ----------
            - file_path (str) : Path to the trajectory file.
            - settings (Settings) : Settings for reading the file.
            
        Returns:
        --------
            - Dict : Dictionary containing the lattice properties.
        """
        # First make sure we have parsed the headers
        if not self.comment_lines:
            self.parse_frames_header(file_path)
            
        # Implementation for LAMMPS format
        # This is a placeholder - actual implementation would parse LAMMPS format
        lattices_dict = {}
        
        # Example implementation (needs to be adapted for actual LAMMPS format)
        for n, line in enumerate(self.comment_lines):
            # Parse lattice information from LAMMPS format
            # This is just a placeholder - actual implementation would be different
            lattices_dict[n] = {
                "lxx": 0.0,
                "lxy": 0.0,
                "lxz": 0.0,
                "lyx": 0.0,
                "lyy": 0.0,
                "lyz": 0.0,
                "lzx": 0.0,
                "lzy": 0.0,
                "lzz": 0.0
            }
            
        return lattices_dict

    def read_atoms(self, file_path: str, settings: Settings) -> Dict:
        """
        Read atomic positions from the trajectory file.
        
        Parameters
        ----------
            - file_path (str) : Path to the trajectory file.
            - settings (Settings) : Settings for reading the file.
            
        Returns:
        --------
            - Dict : Dictionary containing the atomic positions.
        """
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Implementation for LAMMPS format
        # This is a placeholder - actual implementation would parse LAMMPS format
        atoms_dict = {}
        
        # Example implementation (needs to be adapted for actual LAMMPS format)
        # In a real implementation, you would parse the LAMMPS file format
        # and extract atomic positions
        
        return atoms_dict
