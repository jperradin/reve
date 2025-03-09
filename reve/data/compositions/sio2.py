from typing import ClassVar
from .composition import Composition

class SiO2Composition(Composition):
    """SiO2 composition with specific properties for silica systems."""

    name: ClassVar[str] = "SiO2" # Class variable for registration

    def __init__(self) -> None:
        super().__init__(self.name, elements=["Si", "O"])
        self.register_properties()

    def register_properties(self) -> None:
        # Register components for SiO2 system
        self.register_element_component("Si", "NetworkFormer")
        self.register_element_component("O", "NetworkBridge")
        self.register_element_cutoff("Si", "O", 2.30)
        self.register_element_cutoff("O", "O", 3.05)
        self.register_element_cutoff("Si", "Si", 3.50)
        
        # Register system-specific properties
        # Coordination numbers
        self.register_element_property("Si", "coordination", [4, 5, 6, 7, 8])
        self.register_element_property("O", "coordination", [2, 3, 4])
        
        # Average angles
        self.register_element_property("Si", "average_angle", ['OSiO', 'SiSiSi'])
        self.register_element_property("O", "average_angle", ['SiOSi', 'OOO'])
        
        # Average distances
        self.register_element_property("Si", "average_distance", ['SiO', 'SiSi'])
        self.register_element_property("O", "average_distance", ['OO'])
        
        # Pair distribution functions
        self.register_element_property("Si", "pair_distribution_function", ['SiO', 'SiSi'])
        self.register_element_property("O", "pair_distribution_function", ['OO'])
        
        # Bond angle distributions
        self.register_element_property("Si", "bond_angle_distribution", ['OSiO', 'SiSiSi'])
        self.register_element_property("O", "bond_angle_distribution", ['SiOSi', 'OOO'])
        
        # Polyhedricity
        self.register_element_property("Si", "polyhedricity", ['tetrahedricity', 'pentaedricity', 'octahedricity'])
        
        # Number of connectivities
        self.register_element_property("Si", "connectivities", ['number_of_corners', 'number_of_edges', 'number_of_faces'])