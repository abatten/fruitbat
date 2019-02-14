|PyPI| |Python| |Travis| |Docs|

Fruitbat
========

**Fruitbat** is an open source python package used to estimate the redshift of 
Fast Radio Burst (FRB) from its dispersion measure. 

Installation
------------

Run the following to install::

    pip install fruitbat

Or you can clone this repository::
    
    git clone https://github.com/abatten/fruitbat

Usage
-----
To calculate the redshift from a dispersion measure the main function to use
is estimate.redshift_::

    >>> import fruitbat
    >>> fruitbat.estimate.redshift(845)
    2.0

Using estimate.redshift_ you can provide a method and/or a cosmology that you
want to assume.

::

    >>> fruitbat.estimate.redshift(845, method="inoue2004", cosmology="planck2018")
    1.0131760185322125

    >>> fruitbat.estimate.redshift(845, method="zhang2018", cosmology="planck2018+bao")
    1.0131760185322125

.. _estimate.redshift: https://fruitbat.readthedocs.io/en/latest/docstrings/fruitbat.Estimate.html#fruitbat.estimate.redshift


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
