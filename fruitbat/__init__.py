"""
Fruitbat is a package designed for estimating the redshift of FRBs.
"""

from __future__ import absolute_import, print_function, division

from . import utils
from . import cosmologies
from . import plot
from ._frb import Frb
from . import table
from .methods import add_method, reset_methods, avaliable_methods
from .cosmologies import (add_cosmology, reset_cosmologies,
                          avaliable_cosmologies)
from .__version__ import __version__

__name__ = "fruitbat"
__author__ = "Adam Batten (@abatten)"
__all__ = ['Frb', 'methods', 'cosmologies', 'plot', 'table', 'utils',
           'add_method', 'reset_methods', 'avaliable_methods',
           'add_cosmology', 'reset_cosmologies',
           'avaliable_cosmologies']
