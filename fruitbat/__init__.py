"""
Fruitbat
========
Fruitbat is a package designed to assist the estimation of redshifts, energies
and the galactic dispersion measure of fast radio bursts.

"""

from __future__ import absolute_import, print_function, division

from . import utils
from . import methods
from . import cosmologies
from . import table
from . import plot
from ._frb import Frb
from . import catalogue
from .methods import add_method, reset_methods, available_methods
from .cosmologies import (add_cosmology, reset_cosmologies,
                          available_cosmologies)
from .__version__ import __version__

__name__ = "fruitbat"
__author__ = "Adam Batten (@abatten)"
__all__ = ['Frb', 'methods', 'cosmologies', 'plot', 'table', 'utils',
           'catalogue', 'add_method', 'reset_methods', 
           'available_methods', 'add_cosmology', 'reset_cosmologies',
           'available_cosmologies']
