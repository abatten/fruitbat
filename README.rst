|PyPI| |Python| |License| |Travis| |Docs| |CodeCov| |JOSS| |ASCL|

|Logo|

*FRUITBAT* is an open source python package used to estimate the redshift of 
Fast Radio Bursts (FRB) from their dispersion measure. *FRUITBAT* combines 
various dispersion measure (DM) and redshift relations with the YMW16 galactic 
dispersion measure model into a single easy to use API. 

Documentation
-------------
The documentation for *FRUITBAT* can be found at https://fruitbat.readthedocs.io/.

Installation
------------
You can install the latest release of *FRUITBAT* from PyPi_ by running 
the following::

    pip install fruitbat

You can install the latest development version of *FRUITBAT* by cloning 
this repository::
    
    git clone https://github.com/abatten/fruitbat
    cd fruitbat
    pip install .

If you are installing the latest development version of *FRUITBAT* then you 
will also need to install git-lfs. Instructions for installing git-lfs for
your operating system can be found here_.

Linux Users
***********
If you are installing *FRUITBAT* on a linux machine you may see this 'error':
``ERROR: Failed building wheel for pyymw16``. This does not mean the installation
failed. The C++ bindings were compiled using MacOS and needed to be recompiled
for your machine. The installation process does this for you. You should still
be able to run *FRUITBAT* normally.

.. _PyPi: https://pypi.python.org/pypi/fruitbat 
.. _here: https://help.github.com/en/articles/installing-git-large-file-storage


Requirements
------------
Below are the listed requirements for running *FRUITBAT* and the purpose for
each requirement.

 - numpy: Array manipulation

 - scipy: Modules for integration and interpolation

 - astropy: Modules for cosmology, coordinates, constants and units

 - matplotlib: Modules for plotting

 - pandas: Reading csv files from FRBCAT

 - pyymw16: Python wrapper for YMW16 galactic dispersion measure model.

 - e13tools: Utility tools for writing docstrings.

Usage
-----
If you want to get started using *FRUITBAT* there is a `Getting Started`_ 
section of the documentation made just for you! Otherwise the tl;dr is the
following:

Most of the calculations will be centred around the `Frb class`_. You can
can define an instance of the `Frb class`_ with a dispersion measure. 
To calculate the redshift of the FRB use the method 
`calc_redshift`_.

::

    >>> import fruitbat
    >>> FRB121102 = fruitbat.Frb(557, dm_excess=369)
    >>> FRB121102.calc_redshift()
    <Quantity 0.37581945>
    
The `calc_redshift`_ function can also be passed a method and/or a cosmology.
The method will specify which DM-redshift relation to assume and the cosmology
will specify which cosmology to assume.

::

    >>> FRB121102.calc_redshift(method="Zhang2018", cosmology="Planck18")
    <Quantity 0.42166019>

It is also possible to specify the coordinates of the burst and use the 
`calc_dm_galaxy`_ function to calculate the DM contribution from the Milky Way
using the YMW16 or NE2001 galactic electron distribution model. Performing 
`calc_dm_galaxy`_ will automatically calculate the excess dispersion measure 
for the redshift calculation.

::

    >>> FRB190222 = fruitbat.Frb(500, raj="12:34:43.5", decj="2:34:15.2")
    >>> FRB190222.calc_dm_galaxy()
    <Quantity 22.43696785 pc / cm3>
    >>> FRB190222.calc_redshift()
    <Quantity 0.4808557>

.. _Frb class: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Frb.html
.. _calc_redshift: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Frb.html#fruitbat.Frb.calc_redshift
.. _calc_dm_galaxy: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Frb.html#fruitbat.Frb.calc_dm_galaxy
.. _Getting Started: https://fruitbat.readthedocs.io/en/latest/user_guide/getting_started

Issues and Contributing
-----------------------
If there is a feature of *FRUITBAT* that currently does not exist, but you
would like it to, you can contribute by openning a `Github Issue`_ and 
outlining the feature. Similar to contributing, if you find a problem with
*FRUITBAT* or are having difficulties using *FRUITBAT* please do not 
hesitate to open a `Github Issue`_.

.. _Github Issue: https://github.com/abatten/fruitbat/issues

Referencing Fruitbat
--------------------

If you use *FRUITBAT* in your research, we would like it if you could add an 
acknowledgement statement "Some of the results of this paper have been derived
using the *FRUITBAT* package" and reference `our paper`_.

.. _our paper: https://ui.adsabs.harvard.edu/abs/2019JOSS....4.1399B/abstract

::

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




.. |Logo| image:: logo/fruitbat_logo.svg
    :alt: Fruitbat Logo

.. |PyPI| image:: https://img.shields.io/pypi/v/fruitbat.svg?label=PyPI
    :target: https://pypi.python.org/pypi/fruitbat
    :alt: PyPI - Latest Release

.. |Python| image:: https://img.shields.io/pypi/pyversions/fruitbat.svg?label=Python
    :target: https://pypi.python.org/pypi/fruitbat
    :alt: PyPI - Python Versions

.. |Travis| image:: https://travis-ci.com/abatten/fruitbat.svg?branch=master
    :target: https://travis-ci.com/abatten/fruitbat

.. |Docs| image:: https://readthedocs.org/projects/fruitbat/badge/?version=latest
    :target: https://fruitbat.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |CodeCov| image:: https://codecov.io/gh/abatten/fruitbat/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/abatten/fruitbat
    :alt: Code Coverage

.. |License| image:: https://img.shields.io/pypi/l/fruitbat.svg?colorB=purple&label=License
    :target: https://github.com/abatten/fruitbat/raw/master/LICENSE
    :alt: PyPI - License

.. |JOSS| image:: http://joss.theoj.org/papers/634bb69f2445c7457bea5dbc0b83e650/status.svg
    :target: http://joss.theoj.org/papers/634bb69f2445c7457bea5dbc0b83e650
    :alt: JOSS Review Status

.. |ASCL| image:: https://img.shields.io/badge/ascl-1911.010-blue.svg?colorB=262255"
    :target: http://ascl.net/1911.010
    :alt: ascl:1911.010
