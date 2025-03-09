import numpy as np
from typing import List, TextIO, Tuple, Dict, Optional

from..io.base_reader import BaseReader  

class XYZReader(BaseReader):
    """Reader for XYZ trajectory files."""
    def supports_file_format(self, filename: str) -> bool:
        return filename.lower().endswith(".xyz")

    def count_frames(self, filename: str) -> int:
        with open(filename, 'r') as file:
            lines = file.readlines()
        file.close()

        # number of frames = number of comment line
        comments = [line for line in lines if 'Lattice' in line]
        return len(comments)


    def _read_header(self, file: TextIO) -> Tuple[int, Dict]:
        """Reads the XYZ header (num_atoms and comment line)."""
        try:
            num_atoms = int(file.readline())
            comment = file.readline().strip()
            lattice = self._get_lattice(comment)  # Get lattice from comment
            metadata = {"comment": comment, "lattice": lattice, "num_atoms": num_atoms, "num_frames": self.count_frames(file.name)}
            return num_atoms, metadata
        except (ValueError, IOError) as e:  # Catch potential errors
            raise ValueError(f"Invalid XYZ header format: {e}")

    def _read_frame_data(self, file: TextIO, num_atoms: int) -> Optional[List[Dict]]:
        """Reads atom data for a single XYZ frame. Returns None at EOF."""
        frame_data = []
        for _ in range(num_atoms):
            line = file.readline()
            if not line:  # Explicitly check for EOF *before* splitting
                return None  # Signal end of file
            parts = line.split()
            # Check if we have enough parts *after* splitting, but *before* accessing
            if len(parts) < 4:
                return None  # Malformed line (likely EOF or corrupt file)
            try:
                symbol = parts[0]
                x, y, z = map(float, parts[1:4])
                frame_data.append({"symbol": symbol, "position": np.array([x, y, z], dtype=np.float64)})
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping malformed line in XYZ frame data: {e}")
                return None

        return frame_data


    def _read_footer(self, file: TextIO) -> Dict:
        """XYZ files typically don't have footers."""
        return {}

    def _get_lattice(self, comment: str) -> Optional[np.ndarray]:
        """Extracts the lattice matrix from the XYZ comment line, if present."""
        try:
            # Find the Lattice parameter in the comment
            if "Lattice=" not in comment:
                return None  # Or np.eye(3) if you want a default

            # Extract the quoted part
            parts = comment.split('Lattice="')[1].split('"')[0]

            # Split the string into 9 values
            values = list(map(float, parts.split()))

            # Check if we have exactly 9 values for a 3x3 matrix
            if len(values) != 9:
                return None # Or raise an exception

            # Create the 3x3 lattice matrix
            return np.array(values).reshape((3, 3))

        except (IndexError, ValueError):
            return None  # Or raise an exception