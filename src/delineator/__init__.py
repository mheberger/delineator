# src/delineator/__init__.py
import logging
from .core import delineate
from .settings import DelineatorConfig
from .util import write_outputs

__version__ = "0.1.0"

logging.getLogger("delineator").addHandler(logging.NullHandler())
