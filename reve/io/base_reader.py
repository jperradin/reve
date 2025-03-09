from abc import ABC, abstractmethod
from typing import TextIO, List, Tuple, Optional, Dict, Generator
import numpy as np

from ..core.frame import Frame
from ..core.atom import Atom

class BaseReader(ABC):
    """Abstract base class for trajectory readers."""

    @abstractmethod
    def _read_header(self, file: TextIO) -> Tuple[int, Dict]:
        """Reads the header of the trajectory file.

        Args:
            file: File object opened in text mode ('r').

        Returns:
            Tuple: (number of atoms, metadata dictionary).
        """
        pass

    @abstractmethod
    def _read_frame_data(self, file: TextIO, num_atoms: int) -> Optional[List[Dict]]:
        """Reads the data for a single frame.  Crucially, returns None
        when the end of the file is reached.

        Args:
            file: File object opened in text mode ('r').
            num_atoms: Number of atoms to read for this frame.

        Returns:
            Optional[List[Dict]]: List of dictionaries, one for each atom, or None
                if the end of the file has been reached.
        """
        pass


    @abstractmethod
    def _read_footer(self, file: TextIO) -> Dict:
        """Reads the footer of the trajectory file.

        Args:
            file: File object opened in text mode ('r').

        Returns:
            Dict: Metadata dictionary (can be empty).
        """
        pass

    @abstractmethod
    def supports_file_format(self, filename: str) -> bool:
        """Checks if the reader supports the given file format."""
        pass

    def read(self, filename: str) -> Tuple[List[Frame], Dict]:
        """Reads the entire trajectory file, returning all frames at once."""
        frames, metadata = [], {}
        try:
            with open(filename, 'r') as file:
                num_atoms, header_metadata = self._read_header(file)
                metadata.update(header_metadata)

                while True:  # Keep reading until _read_frame_data returns None
                    frame_data = self._read_frame_data(file, num_atoms)
                    if frame_data is None:  # Check for None (end of file)
                        break

                    atoms = [
                        Atom(
                            symbol=atom_data['symbol'],
                            position=atom_data['position'],
                            atom_id=atom_data.get('atom_id'),
                            velocity=atom_data.get('velocity'),
                            mass=atom_data.get('mass'),
                            correlation_length=atom_data.get('correlation_length')
                        )
                        for atom_data in frame_data
                    ]
                    lattice = metadata.get("lattice", np.eye(3))
                    #Crucially, frame_id is determined by the ORDER in the file
                    frame = Frame(atoms=atoms, lattice=lattice, frame_id=len(frames))
                    frames.append(frame)

                footer_metadata = self._read_footer(file)
                metadata.update(footer_metadata)
        except (FileNotFoundError, ValueError, IOError) as e: #Catch IOErrors
            print(f"Error reading file {filename}: {e}")
            return [], {}  # Return empty lists on error
        except Exception as e: #Catch all exceptions.
            print(f"An unexpected error occurred: {e}")
            return [], {}
        return frames, metadata


    def iter_frames(self, filename: str) -> Generator[Frame, None, None]:
      """
      Iterates through the frames of the trajectory file, yielding one Frame at a time.
      This is a generator, so it avoids loading the entire trajectory into memory.
      """
      try:
        with open(filename, 'r') as file:
            num_atoms, header_metadata = self._read_header(file)
            frame_id = 0  # Keep track of the frame ID

            while True:
                frame_data = self._read_frame_data(file, num_atoms)
                if frame_data is None: # Check for None (end-of-file)
                    break

                atoms = [
                    Atom(
                        symbol=atom_data['symbol'],
                        position=atom_data['position'],
                        identifier=atom_data.get('identifier'),
                        velocity=atom_data.get('velocity'),
                        mass=atom_data.get('mass'),
                        correlation_length=atom_data.get('correlation_length')
                    )
                    for atom_data in frame_data
                ]
                lattice = header_metadata.get("lattice", np.eye(3))
                frame = Frame(atoms=atoms, lattice=lattice, frame_id=frame_id)
                yield frame
                frame_id += 1

            # We *could* technically read footer data here, BUT
            # it's generally better to read the footer only once,
            # either with the `read` method or with a separate
            # `read_metadata` method.  Avoid reading the footer
            # repeatedly during iteration.

      except (FileNotFoundError, ValueError, IOError) as e:
        print(f"Error reading file {filename}: {e}")
        return  # Stop iteration on error
      except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return


    def read_metadata(self, filename: str) -> Dict:
      """Reads only the header and footer metadata, without reading frame data."""
      metadata = {}
      try:
        with open(filename, 'r') as file:
            _, header_metadata = self._read_header(file)
            metadata.update(header_metadata)
            # Now, seek to the end and read the footer.  This is a bit tricky
            # because we don't know *exactly* where the footer starts without
            # potentially parsing the entire file.  We'll make an assumption
            # that _read_footer() can handle being called anywhere in the file
            # and will correctly find the footer (if any).  This will need
            # to be handled carefully in the concrete Reader implementations.
            file.seek(0) # reset to the beginning of the file
            footer_metadata = self._read_footer(file)
            metadata.update(footer_metadata)

      except (FileNotFoundError, ValueError, IOError) as e:
            print(f"Error reading metadata from {filename}: {e}")
            return {} #Return empty dict.
      except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return {} # Or raise

      return metadata