import numpy as np
from dataclasses import dataclass, field
from typing import Optional

from ..utils.constants import ATOMIC_MASSES, CORRELATION_LENGTH
from ..utils.geometry import wrap_position

@dataclass(slots=True, order=True)
class Atom:
    """ Reprensation of an atom 
    
    Attributes:
    -----------
    symbol : str
        Symbol of the atom
    atom_id : int
        Id of the atom (unique identifier, autoincremented)
    position : np.ndarray
        Position of the atom
    mass : float
        Mass of the atom (Optional)
    velocity : np.ndarray
        Velocity of the atom (Optional)
    """
    symbol: str
    atom_id: int
    position: np.ndarray = field(compare=False, repr=False)
    velocity: Optional[np.ndarray] = field(default=None, compare=False, repr=False)
    mass: Optional[float] = field(default=None, compare=True, repr=True)
    correlation_length: Optional[float] = field(default=None, compare=True, repr=True)

    _next_id = 1

    def __post_init__(self):
        """ Initialisation after object creation """

        if self.position is None:
            object.__setattr__(self, 'position', np.zeros(3))
        
        if self.velocity is None:
            object.__setattr__(self, 'velocity', np.zeros(3))

        if self.atom_id is None:
            object.__setattr__(self, 'atom_id', Atom._next_id)
            Atom._next_id += 1
        
        if self.mass is None:
            object.__setattr__(self, 'mass', ATOMIC_MASSES.get(self.symbol, 0.0))            

        if self.correlation_length is None:
            object.__setattr__(self, 'correlation_length', CORRELATION_LENGTH.get(self.symbol, 0.0))

    @staticmethod
    def wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
        """ Wrap position in a periodic box defined by the lattice """
        return wrap_position(position, lattice)

    def __str__(self) -> str:
        return f"Atom {self.atom_id} ({self.symbol})"

    def __repr__(self) -> str:
        return self.__str__()
