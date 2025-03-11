"""
REVE - Realistic Environment for Vitreous Exploration
(also verre en verlan)

A package for working with large-scale molecular dynamics trajectories
with a focus on memory efficiency and performance.
"""

# Import core components
from .core.atom import Atom
from .core.frame import Frame
from .core.system import System

# Import IO components
from .io.reader_factory import ReaderFactory

# Import settings
from .settings import Settings

# Import main function
from .main import main

from .utils import print_title

__version__ = "0.1.0"

print_title(__version__)
