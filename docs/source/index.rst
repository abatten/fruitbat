|Logo|

Fruitbat Documentation
----------------------
``Fruitbat`` is an open source Python package designed to assist the estimation of redshifts, 
energies and the galactic dispersion measure contributions of fast radio bursts.

``Fruitbat`` generates and utilises 'look-up' tables of existing dispersion 
measure-redshift relations found in the literature (`Ioka 2003`_, `Inoue 2004`_, 
`Zhang 2018`_) in conjunction with parameters from both the WMAP and Planck 
missions. ``Fruitbat`` also utilises the YMW16 galactic dispersion measure model
to estimate the dispersion measure contribution due to the Milky Way. However it
is also possible to use the NE2001 model if the python port has been installed
(See the `NE2001 installation instructions`_).

As a user you can independantly choose the dispersion measure-reshift relation 
and the cosmological parameters, or define your own relation, create new cosmologies and
generate custom look-up tables for ``fruitbat``.

``Fruitbat`` is installable via ``pip`` (see `Getting Started`_) or the source code is made avaliable here_.

If you use **Fruitbat** in your research, please add the acknowledgement statement
"Some of the results of this paper have been derived using the ``fruitbat`` package"
and cite the JOSS paper.

.. code-block:: none

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

.. _here: https://github.com/abatten/fruitbat
.. _PyPI: https://pypi.org/project/fruitbat
.. _Getting Started: https://fruitbat.readthedocs.io/en/latest/user_guide/getting_started.html
.. _NE2001 installation instructions: https://fruitbat.io/en/latest/user_guide/ne2001_installation.html

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/getting_started
   user_guide/using_fruitbat
   user_guide/method_and_cosmology

.. toctree::
   :maxdepth: 2

   user_guide/faq
   user_guide/ne2001_installation
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
    :target: https://github.com/abatten/fruitbat 

.. _Ioka 2003: https://adsabs.harvard.edu/abs/2003ApJ...598L..79I

.. _Inoue 2004: https://adsabs.harvard.edu/abs/2004MNRAS.348..999I 

.. _Zhang 2018: https://adsabs.harvard.edu/abs/2018ApJ...867L..21Z 

