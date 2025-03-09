
from dataclasses import dataclass
from typing import List, Optional, Dict
import numpy as np
from .element_registry import ElementRegistry


@dataclass(slots=True)
class Atom:
    id: int
    element: str
    position: np.ndarray
    velocity: Optional[np.ndarray] = None
    coordination: Optional[int] = None

    # Call if needed all characteristics of the element from the ElementRegistry.
    @property
    def get_properties(self) -> Dict:
        return ElementRegistry.get_properties(self.element)
    
    @property
    def atomic_number(self) -> int:
        return self.get_properties().get("atomic_number", 0)
    
    @property
    def atomic_mass(self) -> float:
        return self.get_properties().get("atomic_mass", 0.0)
    
    @property
    def correlation_length(self) -> float:
        return self.get_properties().get("correlation_length", 0.0)
    
    @property
    def chemical_symbol(self) -> str:
        return self.get_properties().get("chemical_symbol", "")
    
    @property
    def atomic_name(self) -> str:
        return self.get_properties().get("atomic_name", "")
    
    @property
    def components(self) -> List[str]:
        return self.get_properties().get("components", [])

    def __str__(self):
        return f"Atom({self.id}, {self.element}, {self.position})"
    