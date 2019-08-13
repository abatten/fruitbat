"""
Fruitbat
========
Fruitbat is a package designed to assist the estimation of redshifts, energies
and the galactic dispersion measure of fast radio bursts.

Authors: Adam J. Batten (@abatten)

Short Description
-----------------
This is the ``fruitbat`` package, which allows for easy estimation of FRB
quantities from observed parameters.

A simple example calculate with ``fruitbat`` is to estimate the redshift of 
a FRB from its dispersion measure.

>>> import fruitbat
>>> FRB180525 = fruitbat.FRB(388.1, gl="349", gb="50.7")
>>> dm_milky_way = FRB180525.calc_dm_galaxy(model='ymw16')
>>> frb_redshift = FRB180525.calc_redshift(method='Inoue2004', cosmology="Planck18")
>>> print(frb_redshift, dm_milky_way)
0.3646537633283661 30.56254087351872 pc / cm3

There are more detailed guides for getting started with ``fruitbat`` in the 'Using Fruitbat' section of the online documentation.

https://fruitbat.readthedocs.io/en/latest/user_guide/using_fruitbat.html

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
__author__ = ["Adam Batten (@abatten)", ]
__all__ = ['Frb', 'methods', 'cosmologies', 'plot', 'table', 'utils',
           'catalogue', 'add_method', 'reset_methods', 
           'available_methods', 'add_cosmology', 'reset_cosmologies',
           'available_cosmologies']
