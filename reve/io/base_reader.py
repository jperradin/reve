from abc import ABC, abstractmethod
from typing import Generator, List

from ..core.frame import Frame

class BaseReader(ABC):
    """
    Abstract base class for trajectory readers.  Defines the common interface.
    """
    quiet: bool = False

    def set_quiet(self, quiet: bool) -> None:
        self.quiet = quiet

    @abstractmethod
    def supports_file_format(self, filename: str) -> bool:
        """
        Checks if the reader supports the given file format.

        Args:
            filename: The path to the file.

        Returns:
            bool: True if the reader supports the file format, False otherwise.
        """
        pass

    @abstractmethod
    def read(self, filepath) -> List[Frame]:
        """
        Reads the trajectory data from the given file.

        Args:
            filepath: The path to the trajectory file.

        Returns:
            frames: A list of Frame objects.
        """
        pass

    @abstractmethod
    def get_num_frames(self) -> int:
        """
        Returns the number of frames in the trajectory.

        Returns:
            int: The number of frames.
        """
        pass

    @abstractmethod
    def get_frame(self, frame_index) -> Frame:
        """
        Retrieves a specific frame from the trajectory.

        Args:
            frame_index: The index of the frame to retrieve (0-based).

        Returns:
            Frame: A data structure representing the requested frame.  Should raise an
            IndexError if the frame_index is out of bounds.
        """
        pass

    @abstractmethod
    def iter_frames(self, filename: str, frame_id: int = 0) -> Generator[Frame, None, None]:
        """
        Returns an iterator over the frames in the trajectory.
        """
        pass
    