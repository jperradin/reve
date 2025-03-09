from typing import Dict, List, Any, Optional, ClassVar, Type
from dataclasses import dataclass, field
import importlib
from pathlib import Path
import sys
from ...core.element_registry import ElementRegistry
import numpy as np

@dataclass(slots=True)
class Composition:
    """ Base class for all compositions in the system.

    Uses dataclass with slots=True for memory efficiency when dealing with large systems.
    """
    name: str
    elements: List[str] = field(default_factory=list)
    properties: Dict[str, Any] = field(default_factory=dict)

    # Class variable to store registered compositions   
    _compositions: ClassVar[Dict[str, Type['Composition']]] = {}

    def __init_subclass__(cls, **kwargs) -> None:
        """ Register subclass automatically """
        # Remove the super() call that's causing the error
        if hasattr(cls, "name"):
            Composition._compositions[cls.name] = cls

    @classmethod
    def create(cls, composition_name: str) -> 'Composition':
        """ Factory method to create composition instances

        Args:
            composition_name: Name of the composition to create

        Returns:
            Composition: Composition instance

        Raises:
            ValueError: If the composition name is not registered
        """
        if composition_name in cls._compositions:
            return cls._compositions[composition_name]()

        # If not, try to dynamically import the composition module
        try:
            # Convert composition name to lowercase for file naming convention
            module_name = composition_name.lower()

            # Import the module
            module = importlib.import_module(f"reve.data.compositions.{module_name}")

            # After importing, the element registry will be populated
            # Create a generic composition with registered elements
            elements = ElementRegistry.supported_elements()
            return cls(name=composition_name, elements=elements)

        except ImportError as e:
            raise ValueError(f"Composition '{composition_name}' not found: {e}")

    def register_element_property(self, element: str, property_name: str, property_value: Any) -> None:
        """ Register a property for a specific element in the composition. """
        properties = ElementRegistry.get_properties(element)
        if property_name not in properties:
            properties[property_name] = property_value
        else:
            # Update existing property
            properties[property_name] = property_value

    def register_element_component(self, element: str, component: str) -> None:
        """ Register a component for an element in the composition."""
        ElementRegistry.register_element_component(element, component)
        
    def register_element_cutoff(self, element1: str, element2: str, cutoff: float) -> None:
        """ Register a cutoff distance between two elements in the composition."""
        # Store cutoffs in both elements
        properties1 = ElementRegistry.get_properties(element1)
        properties2 = ElementRegistry.get_properties(element2)

        if "cutoffs" not in properties1:
            properties1["cutoffs"] = {}
        if "cutoffs" not in properties2:
            properties2["cutoffs"] = {}

        properties1["cutoffs"][element2] = cutoff
        properties2["cutoffs"][element1] = cutoff

    def get_elements(self) -> List[str]:
        """ Get the list of elements in this composition """
        return self.elements

    def __str__(self) -> str:
        return f"Composition({self.name}): {', '.join(self.elements)}"