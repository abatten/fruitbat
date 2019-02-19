|PyPI| |Python| |Travis| |Docs|

Fruitbat
========

**Fruitbat** is an open source python package used to estimate the redshift of 
Fast Radio Burst (FRB) from its dispersion measure. 

**NOTE: FRUITBAT IS CURRENTLY IN EARLY DEVELOPMENT AND IS NOT READY FOR SCIENCE YET**

Installation
------------

Run the following to install::

    pip install fruitbat

Or you can clone this repository::
    
    git clone https://github.com/abatten/fruitbat

Usage
-----
To get started you can read the `Getting Started`_ section of the online
documentation.

::

    >>> import fruitbat
    >>> frb = fruitbat.Frb("FRB121102", dm=557, dm_excess=369)
    >>> frb.calc_redshift()
    0.4537827606357084
    


You can provide a method and/or a cosmology that you want to assume.

::

    >>> frb.calc_redshift(method="zhang2018", cosmology="planck2018")
    0.42190276949033323


.. _estimate.redshift: https://fruitbat.readthedocs.io/en/latest/api/fruitbat.Estimate.html#fruitbat.estimate.redshift
.. _Getting Started: https://fruitbat.readthedocs.io/en/latest/user_guide/getting_started

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
