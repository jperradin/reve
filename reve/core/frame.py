import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional

from .atom import Atom
from ..utils.geometry import wrap_positions

@dataclass(slots=True)
class Frame:
    """ 
    Reprensation of a frame of a trajectory

    Attributes:
    -----------
    frame_id : int
        Id of the frame
    atoms : List[Atom]
        List of atoms in the frame
    lattice : Optional[np.ndarray]
        Lattice of the frame
    """
    frame_id: int
    atoms: Optional[List[Atom]] = None
    lattice: Optional[np.ndarray] = None 

    def __post_init__(self):
        """ Initialisation after object creation """
        if self.atoms is not None and not isinstance(self.atoms, list):
            raise TypeError("atoms must be a list of Atoms")
        if self.lattice is not None and not isinstance(self.lattice, np.ndarray):
            raise TypeError("lattice must be a numpy array")
        
    def add_atom(self, atom: Atom) -> None:
        """ Add an Atom object to the list of atoms in the frame """
        if self.atoms is None:
            self.atoms = []
        self.atoms.append(atom)

    def set_lattice(self, lattice: np.ndarray) -> None:
        """ Set the lattice of the frame """
        if lattice.shape != (3, 3):
            raise ValueError("lattice must be a 3x3 numpy array")

        try:
            np.linalg.inv(lattice)
        except np.linalg.LinAlgError:
            raise ValueError("lattice must be a non-singular matrix")

        self.lattice = lattice

    def get_unique_elements(self) -> List[str]:
        """ Get the unique elements in the frame """
        return np.unique([atom.symbol for atom in self.atoms])

    def get_atom_by_id(self, atom_id: int) -> Optional[Atom]:
        """ Get an atom by its id """
        for atom in self.atoms:
            if atom.atom_id == atom_id:
                return atom
        return None

    def get_positions(self) -> np.ndarray:
        """ Get the positions of all atoms in the frame """
        return np.array([atom.position for atom in self.atoms])

    def get_positions_by_element(self) -> Dict[str, np.ndarray]:
        """ Get the positions of all atoms in the frame grouped by element """
        return {atom.symbol: np.array([atom.position for atom in self.atoms if atom.symbol == atom.symbol]) for atom in self.atoms}

    def get_wrapped_positions(self) -> np.ndarray:
        """ Get the wrapped positions of all atoms in the frame """
        return wrap_positions(self.get_positions(), self.lattice)

    def get_wrapped_positions_by_element(self) -> Dict[str, np.ndarray]:
        """ Get the wrapped positions of all atoms in the frame grouped by element """
        return {atom.symbol: wrap_positions(np.array([atom.position for atom in self.atoms if atom.symbol == atom.symbol]), self.lattice) for atom in self.atoms}

    def __len__(self) -> int:
        """ Get the number of atoms in the frame """
        return len(self.atoms)

    def __str__(self) -> str:
        return f"Frame {self.frame_id} (id={self.frame_id})"

    def __repr__(self) -> str:
        return self.__str__()