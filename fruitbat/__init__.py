"""
Fruitbat is a package designed for estimating the redshift of FRBs.
"""

from __future__ import absolute_import, print_function, division

from . import estimate
from . import utils
from . import cosmology
from . import plot
from ._frb import Frb
from .__version__ import __version__

__name__ = "fruitbat"
__author__ = "Adam Batten (@abatten)"
__all__ = ['Frb', 'estimate', 'utils', 'cosmology', 'plot']

