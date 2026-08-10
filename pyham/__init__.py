from .parsers import *
from .ham import *
from .abstractgene import *
from .taxonomy import *
from .genome import *
from .utils import *
from .mapper import *
from .TreeProfile import *
from .iham import *

try:
    from ._version import version as __version__
except ImportError:
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version("pyham")
    except PackageNotFoundError:
        __version__ = "unknown"
