import numpy as np
from typing import List, Optional, Dict, Generator
from tqdm import tqdm

from reve.io.base_reader import BaseReader
from reve.core.frame import Frame
from reve.settings import Settings  # Import the Settings class


class System:
    """
    Manages the atomic system, trajectory data, and interaction with file readers.

    Attributes:
        reader (BaseReader): The file reader used to load data.
        settings (Settings): The settings object containing configuration parameters.
        current_frame (Optional[Frame]): The currently loaded frame.  None if no frame is loaded.
    """

    def __init__(self, reader: BaseReader, settings: Settings):
        """
        Initializes the System object.

        Args:
            reader (BaseReader): The file reader instance to use.
            settings (Settings): The settings object.
        """
        self.reader: BaseReader = reader
        self.settings: Settings = settings
        self.current_frame: Optional[Frame] = None
        self._frame_cache: List[Optional[Frame]] = []  # Cache for faster frame access.
        self._current_frame_index: Optional[int] = None # Index of the current frame
        self._num_frames: Optional[int] = None  # Cache for number of frames

    def load_frame(self, frame_index: int) -> bool:
        """
        Loads a specific frame from the trajectory file.

        Args:
            frame_index (int): The index of the frame to load (0-based).

        Returns:
            bool: True if the frame was successfully loaded, False otherwise.
        """

        if frame_index < 0:
            raise ValueError("Frame index cannot be negative.")

        # Check the range from settings.
        start_frame, end_frame = self.settings.range_of_frames  # Unpack the tuple
        if not (start_frame <= frame_index <= (end_frame if end_frame != -1 else float('inf'))):
            print(f"Frame index {frame_index} is out of range specified in settings ({start_frame}-{end_frame}).")
            return False


        # Check if the frame is in the cache
        if 0 <= frame_index < len(self._frame_cache) and self._frame_cache[frame_index] is not None:
            self.current_frame = self._frame_cache[frame_index]
            self._current_frame_index = frame_index
            return True

        #If not in cache, we use the iterator.
        for frame in self.iter_frames():
            if frame.frame_id == frame_index:
                self.current_frame = frame
                self._current_frame_index = frame_index
                self._add_frame_to_cache(frame) # Add the frame to the cache.
                return True

        print(f"Frame {frame_index} not found in trajectory.")
        return False
    
    def _add_frame_to_cache(self, frame: Frame) -> None:
        """Adds a frame to the cache, managing cache size."""
        if frame.frame_id >= len(self._frame_cache):
            # Extend the cache with None values up to the required index
            self._frame_cache.extend([None] * (frame.frame_id - len(self._frame_cache) + 1))
        self._frame_cache[frame.frame_id] = frame

    def _get_cached_frames(self) -> List[Optional[Frame]]:
        """ Returns the frames in the cache. """
        return self._frame_cache


    def get_frame(self, frame_index: int) -> Optional[Frame]:
        """Retrieves a specific frame, loading it if necessary.

        Args:
            frame_index: The index of the frame to retrieve.

        Returns:
            The Frame object, or None if the frame could not be loaded.
        """
        if self.load_frame(frame_index):
            return self.current_frame
        return None

    def get_num_frames(self) -> int:
        """
        Gets the total number of frames in the trajectory.
        If the value is already calculated, it's directly returned.
        Otherwise, it will be calculated using the method iter_frames.

        Returns:
           int: The total number of frames, 0 if an error occurs
        """
        # First, check if we already calculated the number of frames
        return self.reader.get_num_frames(self.settings.file_location)

    def iter_frames(self) -> Generator[Frame, None, None]:
        """
        Iterates through the frames of the trajectory, yielding one Frame at a time.
        This is a generator, avoiding loading the entire trajectory into memory.
        It respects the range of frames defined in settings.

        Yields:
            Frame: The next Frame object in the trajectory.
        """
        start_frame, end_frame = self.settings.range_of_frames
        # Use the reader's iterator to avoid loading everything into memory.
        frame_generator = self.reader.iter_frames(self.settings.file_location)

        # Skip initial frames if start_frame > 0
        for _ in range(start_frame):
            try:
                next(frame_generator) # Consume frames, and increment frame_id
            except StopIteration:
                return # Stop if the file has fewer frames than start_frame
            

        # Yield frames within the specified range
        counter = 0
        for frame_id, frame in enumerate(frame_generator, start=start_frame):
            if end_frame != -1 and frame_id > end_frame:
                break  # Stop if we've reached the end frame
            yield frame
            counter += 1


    def __iter__(self) -> 'System':
        """
        Make the System object itself iterable.
        """
        # Reset the frame index for iteration.
        self._current_frame_index = self.settings.range_of_frames[0]
        return self


    def __next__(self) -> Frame:
        """
        Returns the next frame during iteration.
        """

        if self._current_frame_index is None:  # First call to next
             self._current_frame_index = self.settings.range_of_frames[0] # Initialize if needed.

        start_frame, end_frame = self.settings.range_of_frames

        if end_frame != -1 and self._current_frame_index > end_frame: #check if end frame is reach
            raise StopIteration

        if self._current_frame_index < self.get_num_frames():
            if self.load_frame(self._current_frame_index):
                self._current_frame_index += 1
                return self.current_frame  # type: ignore
            else: # If load frame return False
                raise StopIteration
        else: # if current frame index is greater than the number of frames.
            raise StopIteration