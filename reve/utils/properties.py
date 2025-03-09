from dataclasses import dataclass, field
from typing import Dict

@dataclass(frozen=True)
class AtomProperties:
    symbol: str
    properties : Dict

ATOM_PROPERTIES_DATA = [
    AtomProperties('Si', {}),
    AtomProperties('O', {}),
    AtomProperties('H', {}),
    AtomProperties('Na', {}),
    AtomProperties('B', {}),
    # TODO: implement additional properties when needed
]

__all__ = [
    'AtomProperties',
    'ATOM_PROPERTIES_DATA'
]