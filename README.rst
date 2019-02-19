|PyPI| |Python| |Travis| |Docs|

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

Usage
-----
If you want to get started using **fruitbat** there is a `Getting Started`_ 
section of the documentation made just for you! Otherwise the tl;dr is the
following:

Most of the calculations will be centred around the `Frb class`_. You can
can define an instance of the `~Frb` with a name and a dispersion 
measure. To calculate the redshift of the FRB use the method 
`calc_redshift_` ::

    >>> import fruitbat
    >>> frb = fruitbat.Frb("FRB121102", dm=557, dm_excess=369)
    >>> frb.calc_redshift()
    0.4537827606357084
    
You can provide a method and/or a cosmology that you want to assume.

::

    >>> frb.calc_redshift(method="zhang2018", cosmology="planck2018")
    0.42190276949033323


.. _Frb class: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Frb.html
.. _calc_redshift: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Frb.html#fruitbat.Frb.calc_redshift
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
