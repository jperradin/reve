
from typing import Generator
import numpy as np

from .base_reader import BaseReader
from ...core.frame import Frame
from ...core.atom import Atom

class XYZReader(BaseReader):
    """
    Abstract base class for trajectory readers.  Defines the common interface.
    """
    verbose: bool = True

    def supports_file_format(self, filename: str) -> bool:
        return filename.lower().endswith('.xyz')

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
                if 'Lattice' in line:
                    num_frames += 1
        f.close()
        return num_frames

    def get_frame(self, frame_index) -> Frame:
        """
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
                    27216
                    Lattice="65.0 0.0 0.0 0.0 65.0 0.0 0.0 0.0 65.0" Properties=...
                    """
                    try:
                        num_atoms = int(f.readline())
                    except ValueError:
                        # End of file or invalid line.  Stop the generator.
                        break  # Important: Exit the while loop

                    comment = f.readline().strip()  # Read and discard the comment line.
                    # Extract lattice in header
                    lattice = comment.split('\"')[1]
                    l = np.array(lattice.split(), dtype=float)
                    lx = np.array([l[0],l[1], l[2]])
                    ly = np.array([l[3],l[4], l[5]])
                    lz = np.array([l[6],l[7], l[8]])
                    lattice = np.array([lx, ly, lz])

                    atoms = []
                    for i in range(num_atoms):
                        # Read lines in file
                        line = f.readline().strip()
                        if not line:  # Check for unexpected EOF within a frame
                            raise ValueError("Unexpected end of file within a frame.")
                        parts = line.split()
                        
                        # Create Atom
                        symbol, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
                        atom = Atom(atom_id=i, symbol=symbol, position=np.array([x, y, z]))
                        if len(parts) > 4:
                            atom.other = parts[4:]
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