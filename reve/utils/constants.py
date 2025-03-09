from dataclasses import dataclass, field
from typing import Dict

@dataclass(frozen=True)
class AtomicMass:
    symbol: str
    mass: float

@dataclass(frozen=True)
class CorrelationLength:
    symbol: str
    length: float

ATOMIC_MASSES_DATA = [
    AtomicMass('Si', 28.0855),
    AtomicMass('O', 15.9994),
    AtomicMass('H', 1.0080),
    AtomicMass('Na', 22.9897),
    AtomicMass('B', 10.8115),
    # TODO: implement more elements when needed
]

CORRELATION_LENGTH_DATA = [
    CorrelationLength('Si', 4.149),
    CorrelationLength('O', 5.803),
    CorrelationLength('H', -3.7390),
    CorrelationLength('Na', 3.63),
    CorrelationLength('B', 5.30),
    # TODO: implement more elements when needed
]

ATOMIC_MASSES: Dict[str, float] = {am.symbol: am.mass for am in ATOMIC_MASSES_DATA}
CORRELATION_LENGTH: Dict[str, float] = {cl.symbol: cl.length for cl in CORRELATION_LENGTH_DATA}

def get_atomic_mass(symbol: str) -> float:
    return ATOMIC_MASSES.get(symbol, 0.0) 

def get_correlation_length(symbol: str) -> float:
    return CORRELATION_LENGTH.get(symbol, 0.0) 

__all__ = [
    'ATOMIC_MASSES',
    'CORRELATION_LENGTH',
    'get_atomic_mass',
    'get_correlation_length'
]