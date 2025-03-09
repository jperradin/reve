from typing import List, TextIO, Tuple, Dict, Optional

from ..core.frame import Frame
from .base_reader import BaseReader


class LAMMPSReader(BaseReader):
    """
    Reader for LAMMPS trajectory files
    """

    def _supports_file_format(self, file: TextIO) -> bool:
        """
        Check if the file format is supported by the reader
        """
        return file.name.endswith('.data')

    def _read_header(self, file: TextIO) -> Tuple[int, Dict]:
        """
        Read the header of the trajectory file

        Parameters
        ----------
        file : TextIO
            File object to read from

        Returns
        -------
        Tuple[int, Dict]
            Tuple containing the frame id and the header lines as a metadata dictionary.
        """
        return 0, {}

    def _read_footer(self, file: TextIO) -> Tuple[int, Dict]:
        """
        Read the footer of the trajectory file

        Parameters
        ----------
        file : TextIO
            File object to read from

        Returns
        -------
        Tuple[int, Dict]
            Tuple containing the frame id and the footer lines as a metadata dictionary.
        """
        return 0, {}

    def _read_frame(self, file: TextIO, n_atoms: int) -> Frame:
        """
        Read a frame from the trajectory file

        Parameters
        ----------
        file : TextIO
            File object to read from
        n_atoms : int
            Number of atoms in the frame

        Returns
        -------
        List[Atom]
            List of Atoms in the frame
        """
        return Frame(0, [], np.zeros((3, 3)))