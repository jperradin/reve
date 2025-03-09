from abc import ABC, abstractmethod
from typing import Dict, List, Optional, List, Any, Type, Callable
import numpy as np
from ..core.frame import Frame
from ..settings.settings import Settings


class Analyzer(ABC):
    """
    Base class for analyzers with factory method pattern.
    """

    # Registry to store analyzer types
    _analyzers: Dict[str, Type["Analyzer"]] = {}

    def __init__(self, settings: Settings) -> None:
        self.settings = settings
        self.results = {}

    @classmethod
    def register(cls, analyzer_type: str) -> Callable:
        """ Decorator to register analyzer types. """
        def decorator(analyzer_class: Type["Analyzer"]):
            cls._analyzers[analyzer_type] = analyzer_class
            return analyzer_class
        return decorator

    @classmethod   
    def create(cls, analyzer_type: str, settings: Settings) -> "Analyzer":
        """ Factory method to create analyzer instances. """
        if analyzer_type not in cls._analyzers:
            raise ValueError(f"Unknown analyzer type: {analyzer_type}")
        return cls._analyzers[analyzer_type](settings)
        
    @classmethod
    def available_analyzers(cls) -> List[str]:
        """ Return list of available analyzer types. """
        return list(cls._analyzers.keys())

    @abstractmethod
    def analyze(self, frames: List[Frame]) -> Dict[str, Any]:
        """ Analyze the frames and return results. """
        pass

    @abstractmethod
    def finalize(self) -> Dict[str, Any]:
        """ Finalize the analysis after all frames are processed. """
        pass

    def save_results(self, output_path: Optional[str] = None) -> None:
        """ Save the results to a file. """
        if output_path is None:
            output_path = self.settings._output_directory

        # implementation for saving results
        # implement writer factory in io.writer
        pass