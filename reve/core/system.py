from typing import Dict, Optional, Callable
from .frame import Frame
from ..settings.settings import Settings
from ..data.compositions.composition import Composition

class System:
    """
    Represents a molecular system with multiple frames.
    
    Attributes:
        name (str): Name of the system
        current_frame (Frame): Current active frame
    """
    
    def __init__(self, name: str, settings: Settings) -> None:
        self.name: str = name
        self.current_frame: Frame = None
        self.analyzers: Dict[str, Analyzer] = {}

        # Load composition from settings using the factory
        composition_name = settings.composition.value if hasattr(settings, 'composition') else 'SiO2'
        self.composition = Composition.create(composition_name)
        
    def get_current_frame(self) -> Optional[Frame]:
        """Get the current active frame"""
        return self.current_frame
    
    def set_current_frame(self, frame: Frame) -> None:
        """Set the current active frame"""
        self.current_frame = frame
    
    def __str__(self) -> str:
        return f"System {self.name}:, current frame: {self.current_frame.id}, composition: {self.composition}, analyzers: {self.analyzers}"