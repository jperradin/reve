from typing import List, Generator
from collections import namedtuple
import numpy as np
import os

from .base_reader import BaseReader
from ...core.frame import Frame
from ...core.atom import Atom

FrameIndex = namedtuple('FrameIndex', ['frame_id', 'num_atoms', 'lattice', 'byte_offset'])

class XYZReader(BaseReader):
    """
    Reader for XYZ trajectory files.
    """

    def detect(self, filepath: str) -> bool:
        """
        Detects if the file is supported by this reader.

        Returns:
            bool: True if the file is supported, False otherwise.
        """
        return filepath.lower().endswith('.xyz')

    def scan(self) -> List[FrameIndex]:
        """
        Scans the trajectory file.
        Initializes Frame objects with the chunk locations of each frame.
        Parse the header to store the number of atoms, the lattice and other informations.

        Returns:
            List[FrameIndex]: A list of FrameIndex objects.
        """
        if not self.mmaped_file:
            with open(self.filename, 'rb') as f:
                self.mmaped_file = memoryview(f.read())
        self.frame_offsets = []
        self.frame_sizes = []
        self.num_frames = 0
        self.frame_indices = []

        try:
            file_size = os.path.getsize(self.filename)
            with open(self.filename, 'r') as f:
                offset = 0
                while offset < file_size:
                    # Store position at the start of the frame
                    frame_start = offset

                    # Read number of atoms line
                    num_atoms_line = f.readline().strip()
                    if not num_atoms_line:
                        break # End of file
                    try:
                        num_atoms = int(num_atoms_line)
                    except ValueError:
                        raise ValueError("Number of atoms must be an integer")
                    offset += 1

                    # Read lattice line
                    lattice_line = f.readline().strip()
                    try:
                        lattice_line = lattice_line.split('\"')[1]
                        lattice_line = lattice_line.split()
                        lxx, lxy, lxz = float(lattice_line[0]), float(lattice_line[1]), float(lattice_line[2])
                        lyx, lyy, lyz = float(lattice_line[3]), float(lattice_line[4]), float(lattice_line[5])
                        lzx, lzy, lzz = float(lattice_line[6]), float(lattice_line[7]), float(lattice_line[8])
                        lattice = np.array([[lxx, lxy, lxz], [lyx, lyy, lyz], [lzx, lzy, lzz]])
                    except ValueError:
                        raise ValueError("Lattice must be a 3x3 matrix")
                    offset += 1

                    # Read atom lines
                    for _ in range(num_atoms):
                        l = f.readline()
                        offset += 1
                    
                    # Store position at the end of the frame
                    self.frame_indices.append(FrameIndex(frame_id=self.num_frames, num_atoms=num_atoms, lattice=lattice, byte_offset=frame_start))
                    self.num_frames += 1
                    self.frame_sizes.append(offset - frame_start)
                    self.frame_offsets.append(frame_start)

        except Exception as e:
            raise Exception(f"Error scanning trajectory file: {str(e)}")

        if self.verbose:
            print(f"Scanned {self.num_frames} frames in {self.filename}")

        self.is_indexed = True
        return self.frame_indices

    def parse(self, frame_id: int) -> Generator[Frame, None, None]:
        """
        Parses the trajectory file, get atom data and yields frames.

        Yields:
            Frame: A data structure representing a frame.
        """
        
        frame_index = self.frame_indices[frame_id]
        with open(self.filename, 'r') as f:
            self.seek_to_line(frame_index.byte_offset) 
            num_atoms = frame_index.num_atoms
            lattice = frame_index.lattice

            # Skip 2 header lines
            f.readline()
            f.readline()

            symbols = []
            positions = []
            # Read atom lines
            for _ in range(num_atoms):
                atom_line = f.readline().strip()
                try:
                    parts = atom_line.split()
                    symbol = parts[0]
                    x, y, z = map(float, parts[1:])
                    symbols.append(symbol)
                    positions.append(np.array([x, y, z]))
                except ValueError:
                    raise ValueError("Atom line must have 4 values: symbol, x, y, z")

            data = {
                'symbol': symbols,
                'position': positions
            }
            
            yield Frame(frame_id=frame_id, _data=data, lattice=lattice, atoms=[])
                
                

