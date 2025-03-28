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
from .io.reader.reader_factory import ReaderFactory
# from .io.writer.writer_factory import WriterFactory # TODO: implement writers

# Import settings
from .config.settings import Settings, SettingsBuilder, PDFSettings, BADSettings, AnalysisSettings

# Import main function
from .main import main

from .utils import print_title

__version__ = "0.0.1"

print_title(__version__)
