|Logo|

Fruitbat Documentation
----------------------
**Fruitbat** is an open source Python package designed to assist the estimation of redshifts, 
energies and the galactic dispersion measure contributions of fast radio bursts.

**Fruitbat** generates and utilises 'look-up' tables of existing dispersion measure-redshift relations
found in the literature (`Ioka 2003`_, `Inoue 2004`_, `Zhang 2018`_) in conjunction with parameters 
from both the WMAP and Planck missions. **Fruitbat** also utilised the YMW16 galactic dispersion measure model
to estimate the dispersion measure contribution due to the Milky Way.

As a user you can independantly choose the dispersion measure-reshift relation 
and the cosmological parameters, or define your own relation, create new cosmologies and
generate custom look-up tables for **fruitbat**.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/getting_started
   user_guide/using_fruitbat
   user_guide/method_and_cosmology

.. toctree::
   :maxdepth: 2

   user_guide/faq
   user_guide/guidelines


.. toctree::
    :maxdepth: 2
    :caption: API Documentation

    api/fruitbat.Frb
    api/fruitbat.methods
    api/fruitbat.cosmologies
    api/fruitbat.catalogue
    api/fruitbat.table
    api/fruitbat.plot
    api/fruitbat.utils




.. |Logo| image:: ../../logo/fruitbat_logo.svg 

.. _Ioka 2003: https://adsabs.harvard.edu/abs/2003ApJ...598L..79I

.. _Inoue 2004: https://adsabs.harvard.edu/abs/2004MNRAS.348..999I 

.. _Zhang 2018: https://adsabs.harvard.edu/abs/2018ApJ...867L..21Z 

