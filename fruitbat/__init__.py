"""
Fruitbat
========
Fruitbat is a package designed to assist the estimation of redshifts, energies
and the galactic dispersion measure of fast radio bursts.

Authors: Adam J. Batten (@abatten)

Short Description
-----------------
This is the ``FRUITBAT`` package, which allows for easy estimation of FRB
quantities from observed parameters.

A simple example calculate with *FRUITBAT* is to estimate the redshift of
a FRB from its dispersion measure.

>>> import fruitbat
>>> FRB180525 = fruitbat.FRB(388.1, gl="349", gb="50.7")
>>> dm_milky_way = FRB180525.calc_dm_galaxy(model='ymw16')
>>> frb_redshift = FRB180525.calc_redshift(method='Inoue2004', cosmology="Planck18")
>>> print(frb_redshift, dm_milky_way)
0.3646537633283661 30.56254087351872 pc / cm3

There are more detailed guides for getting started with *FRUITBAT* in the 'Using Fruitbat'
section of the online documentation.

https://fruitbat.readthedocs.io/en/latest/user_guide/using_fruitbat.html


If you use ``FRUITBAT`` in your research, we would like it if you could reference our paper.
You can get the bibtex using the following function:

>>> fruitbat.get_bibtex()

"""

from __future__ import absolute_import, print_function, division

from textwrap import dedent

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

__pkgname__ = "FRUITBAT"
__author__ = ["Adam Batten (@abatten)", ]
__all__ = ['Frb', 'methods', 'cosmologies', 'plot', 'table', 'utils',
           'catalogue', 'add_method', 'reset_methods',
           'available_methods', 'add_cosmology', 'reset_cosmologies',
           'available_cosmologies']




def get_bibtex():
    # Thanks to Ellert van der Velden (@1313e) for this code that was blatently copied from PRISM.
    """
    Prints a string that gives the BibTeX entry for citing the ``FRUITBAT`` paper
    (Batten 2019, JOSS, vol. 4, issue 17, id. 1399).

    """

    # Create string with BibTeX entry
    bibtex = dedent(
        r"""
        @ARTICLE{2019JOSS....4.1399B,
               author = {{Batten}, Adam},
                title = "{Fruitbat: A Python Package for Estimating Redshifts of Fast Radio Bursts}",
              journal = {The Journal of Open Source Software},
             keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - High Energy Astrophysical Phenomena},
                 year = "2019",
                month = "May",
               volume = {4},
               number = {37},
                pages = {1399},
                  doi = {10.21105/joss.01399},
        archivePrefix = {arXiv},
               eprint = {1905.04294},
         primaryClass = {astro-ph.IM},
               adsurl = {https://ui.adsabs.harvard.edu/abs/2019JOSS....4.1399B},
              adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }
        """)

    # Print the string
    print(bibtex.strip())
    return bibtex.strip()

def __cite__():
    return get_bibtex()