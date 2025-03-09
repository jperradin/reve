"""
REVE - Realistic Environment for Vitreous Exploration
(also verre en verlan)

A package for working with large-scale molecular dynamics trajectories
with a focus on memory efficiency and performance.
"""

# Import core components
from .core.atom import Atom
from .core.frame import Frame

# Import IO components
from .io.reader_factory import ReaderFactory

# Import settings
from .settings import Settings

# Import main function
from .main import main

__version__ = "0.1.0"

TITLE = fr'''
                                                 
`7MM"""Mq.  `7MM"""YMM `7MMF'   `7MF'`7MM"""YMM  
  MM   `MM.   MM    `7   `MA     ,V    MM    `7  
  MM   ,M9    MM   d      VM:   ,V     MM   d    
  MMmmdM9     MMmmMM       MM.  M'     MMmmMM    
  MM  YM.     MM   Y  ,    `MM A'      MM   Y  , 
  MM   `Mb.   MM     ,M     :MM;       MM     ,M 
.JMML. .JMM..JMMmmmmMMM      VF      .JMMmmmmMMM 

'''
print(TITLE)
print(f"version \u279c\t {__version__}\n")