import numpy as np
from numba import jit

@jit(nopython=True, cache=True)
def wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """ Wrap position in a periodic box defined by the lattice """
    position = np.ascontiguousarray(position)
    lattice = np.ascontiguousarray(lattice)

    lattice_inv = np.linalg.inv(lattice)
    lattice_inv = np.ascontiguousarray(lattice_inv)

    fractional_position = np.dot(position, lattice_inv)
    fractional_position = np.ascontiguousarray(fractional_position)

    fractional_position -= np.floor(fractional_position)
    fractional_position = np.ascontiguousarray(fractional_position)

    wrapped_position = np.dot(fractional_position, lattice)
    wrapped_position = np.ascontiguousarray(wrapped_position)

    return wrapped_position

@jit(nopython=True, cache=True)
def wrap_positions(positions: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """ Wrap positions in a periodic box defined by the lattice """
    wrapped_positions = np.zeros_like(positions)
    for i in range(positions.shape[0]):
        position = np.ascontiguousarray(positions[i])
        lattice = np.ascontiguousarray(lattice)

        lattice_inv = np.linalg.inv(lattice)
        lattice_inv = np.ascontiguousarray(lattice_inv)

        fractional_position = np.dot(position, lattice_inv)
        fractional_position = np.ascontiguousarray(fractional_position)

        fractional_position -= np.floor(fractional_position)
        fractional_position = np.ascontiguousarray(fractional_position)

        wrapped_position = np.dot(fractional_position, lattice)
        wrapped_position = np.ascontiguousarray(wrapped_position)

        wrapped_positions[i] = wrapped_position

    return wrapped_positions

__all__ = [
    'wrap_position',
    'wrap_positions'
]