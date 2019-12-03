# Ensure compatibility with Python 2
from __future__ import absolute_import, division, print_function, unicode_literals

from .version import __version__
__all__ = ['__version__']

from .interface import solve
__all__ += ['solve']