
from typing import Generator, Optional
import numpy as np
import re

from .base_reader import BaseReader
from ..core.frame import Frame
from ..core.atom import Atom

class LAMMPSReader(BaseReader):
    """
    Abstract base class for trajectory readers.  Defines the common interface.
    """
    quiet: bool = False

    def supports_file_format(self, filename: str) -> bool:
        return filename.lower().endswith('.lammpstrj')

    def read(self, filepath):
        """
        Reads the trajectory data from the given file.

        Args:
            filepath: The path to the trajectory file.

        Returns:
            frames: A list of Frame objects.
        """
        return list(self.iter_frames(filepath))
        

    def get_num_frames(self, filepath: str) -> int:
        """
        Returns the number of frames in the trajectory.

        Returns:
            int: The number of frames.
        """
        num_frames = 0
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'ITEM: ATOMS' in line:
                    num_frames += 1
        f.close()
        return num_frames

    def get_frame(self, frame_index) -> Frame:
        """
        Retrieves a specific frame from the trajectory.
        Retrieves a specific frame from the trajectory.

        Args:
            frame_index: The index of the frame to retrieve (0-based).

        Returns:
            Frame: A data structure representing the requested frame.  Should raise an
            IndexError if the frame_index is out of bounds.
        """
        return list(self.iter_frames(filepath))[frame_index]

    def iter_frames(self, filename: str, frame_id: int = 0) -> Generator[Frame, None, None]:
        """
        Returns an iterator over the frames in the trajectory.
        """
        try:
            with open(filename, 'r') as f:
                while True:  # Keep reading until the file is exhausted.
                    """ Header is like :
                    ITEM: TIMESTEP
                    0
                    ITEM: NUMBER OF ATOMS
                    27216
                    ITEM: BOX BOUNDS pp pp pp
                    2.8250412707498995e+00 7.1174958729248232e+01
                    2.8250412707498995e+00 7.1174958729248232e+01
                    2.8250412707498995e+00 7.1174958729248232e+01
                    ITEM: ATOMS id type x y z vx vy vz c_pot c_peratom[1] c_peratom[2] c_peratom[3] c_peratom[4] c_peratom[5] c_peratom[6]
                    """
                    for skip in range(3):
                        f.readline()
                    # Read number of atoms
                    try:
                        num_atoms = int(f.readline())
                    except ValueError:
                        # End of file or invalid line.  Stop the generator.
                        break  # Important: Exit the while loop

                    # skip nest line
                    f.readline()

                    # Extract lattice in header
                    lxlo, lxhi = f.readline().strip().split()
                    lylo, lyhi = f.readline().strip().split()
                    lzlo, lzhi = f.readline().strip().split()
                    lxx = np.float64(lxhi) - np.float64(lxlo)
                    lyy = np.float64(lyhi) - np.float64(lylo)
                    lzz = np.float64(lzhi) - np.float64(lzlo)

                    lx = np.array([lxx, 0, 0])
                    ly = np.array([0, lyy, 0])
                    lz = np.array([0, 0, lzz])

                    lattice = np.array([lx, ly, lz])

                    comment = f.readline().strip()  # Read and discard the comment line.
                    
                    # Decode comment line by fetching keywords and columns
                    comment = comment.split()[1:]
                    # Test if ATOMS is present
                    if 'ATOMS' not in comment:
                        raise ValueError("Invalid LAMMPS file format: ATOMS keyword not found. Check the file format.")

                    # Remove 'Item:' and 'ATOMS' with -2
                    comment = [c for c in comment if c not in ['Item:', 'ATOMS']]

                    # Mandatory columns if not present, raise error
                    try:
                        type_column = comment.index('type')
                        x_column = comment.index('x')
                        y_column = comment.index('y')
                        z_column = comment.index('z')
                    except ValueError:
                        raise ValueError("Invalid LAMMPS file format: ATOMS keyword not found. Check the file format.")
                    
                    # Optional columns
                    id_column = comment.index('id') if 'id' in comment else None
                    vx_column = comment.index('vx') if 'vx' in comment else None
                    vy_column = comment.index('vy') if 'vy' in comment else None
                    vz_column = comment.index('vz') if 'vz' in comment else None
                    
                    # Extract atoms
                    atoms = []
                    for i in range(num_atoms):
                        line = f.readline().strip()
                        if not line:  # Check for unexpected EOF within a frame
                            raise ValueError("Unexpected end of file within a frame.")
                        parts = line.split()
                        symbol, x, y, z = parts[type_column], float(parts[x_column]), float(parts[y_column]), float(parts[z_column])
                        
                        # Check optional columns
                        if id_column is not None:
                            id = int(parts[id_column])
                        else:
                            id = i
                        if vx_column is not None:
                            vx = float(parts[vx_column])
                        else:
                            vx = 0.0
                        if vy_column is not None:
                            vy = float(parts[vy_column])
                        else:
                            vy = 0.0
                        if vz_column is not None:
                            vz = float(parts[vz_column])
                        else:
                            vz = 0.0
                        
                        # Create atom
                        atom = Atom(atom_id=id, symbol=symbol, position=np.array([x, y, z]), velocity=np.array([vx, vy, vz]))
                        atoms.append(atom)

                    # Create frame and yield it
                    frame = Frame(frame_id=frame_id, atoms=atoms, lattice=lattice)
                    yield frame
                    frame_id += 1

        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.")
            # You might choose to raise the exception here or return an empty generator.
            # Raising is generally preferred, as it clearly signals an error.
            return  # or `raise`

        except ValueError as e:
            print(f"Error: Invalid XYZ file format: {e}")
            return # or `raise`