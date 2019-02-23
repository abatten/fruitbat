|PyPI| |Python| |Travis| |Docs| |CodeCov|

Fruitbat
========

**Fruitbat** is an open source python package used to estimate the redshift of 
Fast Radio Burst (FRB) from its dispersion measure. 

**NOTE: FRUITBAT IS CURRENTLY IN EARLY DEVELOPMENT AND IS NOT READY FOR SCIENCE YET**


Documentation
-------------

The documentation for **fruitbat** can be found at 
https://fruitbat.readthedocs.io/en/latest/index.html.

Installation
------------

Run the following to install::

    pip install fruitbat

Or you can clone this repository::
    
    git clone https://github.com/abatten/fruitbat
    cd fruitbat
    pip install .

Requirements
------------
 - numpy

 - scipy

 - astropy

 - pyymw16

Usage
-----
If you want to get started using **fruitbat** there is a `Getting Started`_ 
section of the documentation made just for you! Otherwise the tl;dr is the
following:

Most of the calculations will be centred around the `Frb class`_. You can
can define an instance of the `Frb class`_ with a name and a dispersion 
measure. To calculate the redshift of the FRB use the method 
`calc_redshift_` ::

    >>> import fruitbat
    >>> frb = fruitbat.Frb("FRB121102", dm=557, dm_excess=369)
    >>> frb.calc_redshift()
    0.4537827606357084
    
You can also provide a method and/or a cosmology that you want to assume.

::

    >>> frb.calc_redshift(method="zhang2018", cosmology="planck2018")
    0.42190276949033323

It is also possible to specify the coordinates of the burst and use the 
`calc_dm_galaxy_` function to calculate the DM contribution from the Milky Way
using the YMW16 model. Performing `calc_dm_galaxy_` will automatically
calculate the excess dispersion measure for the redshift calculation.

::
    >>> FRB190222 = fruitbat.Frb("FRB190222", dm=500, raj="12:34:43.5", decj="2:34:15.2")
    >>> frb.calc_dm_galaxy()
    22.436967849731445
    >>> frb.calc_redshift()
    0.48112390552750095



.. _Frb class: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Frb.html
.. _calc_redshift: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Frb.html#fruitbat.Frb.calc_redshift
.. _calc_dm_galaxy: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Frb.html#fruitbat.Frb.calc_dm_galaxy
.. _Getting Started: https://fruitbat.readthedocs.io/en/latest/user_guide/getting_started


Referencing Fruitbat
--------------------
If you use ``fruitbat`` in your research, we would like it if you could
reference us. :)


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
