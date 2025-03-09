from typing import List, Dict
from dataclasses import dataclass
import numpy as np
from .atom import Atom
from numba import njit, float64, int64 

@dataclass
class Frame:
    id: int
    # Lattice parameters
    lattice: np.ndarray

    # Store bulk data in arrays for efficient operations
    positions: np.ndarray   
    elements: np.ndarray
    atoms_ids: np.ndarray

    _atoms: Dict[int, Atom] = None

    def get_atoms(self) -> Dict[int, Atom]:
        if self._atoms is None:
            self._atoms = {i: Atom(i, self.elements[i], self.positions[i]) for i in range(len(self.positions))}
        return self._atoms

    # Wrap positions at initialization
    @property
    def positions(self) -> np.ndarray:
        return self._positions

    @positions.setter
    def positions(self, value: np.ndarray) -> None:
        positions = np.ascontiguousarray(value)
        lattice = np.ascontiguousarray(self.lattice)
        self._positions = self.wrap_positions(positions, lattice)

    @property
    def elements(self) -> np.ndarray:
        return self._elements

    @elements.setter
    def elements(self, value: np.ndarray) -> None:
        self._elements = value

    @property
    def atoms_ids(self) -> np.ndarray:
        return self._atoms_ids

    @atoms_ids.setter
    def atoms_ids(self, value: np.ndarray) -> None:
        self._atoms_ids = value

    def __str__(self) -> str:
        return f"Frame({self.id}, Lattice {self.lattice}, {len(self.positions)} atoms)"

    def __repr__(self) -> str:
        return self.__str__()

    @staticmethod
    @njit(float64[:, :](float64[:, :], float64[:, :]))
    def wrap_positions(positions: np.ndarray, lattice: np.ndarray) -> np.ndarray:
        # Compute the inverse of the lattice matrix
        lattice_inv = np.linalg.inv(lattice)
        lattice_inv = np.ascontiguousarray(lattice_inv)
        positions = np.ascontiguousarray(positions)
        # Convert positions to fractional coordinates
        s_positions = np.dot(positions, lattice_inv)
        s_positions = np.ascontiguousarray(s_positions)

        # Wrap positions inside the unit cell
        s_positions -= np.floor(s_positions)
        s_positions = np.ascontiguousarray(s_positions)
        lattice = np.ascontiguousarray(lattice)
        
        # Transform positions back to the Cartesian coordinate system
        wrapped_positions = np.dot(s_positions, lattice)
        wrapped_positions = np.ascontiguousarray(wrapped_positions)
        
        return wrapped_positions