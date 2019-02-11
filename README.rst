|Travis| |Docs|

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

The main function you want to use is `dm_to_redshift()`::

    >>> import fruitbat
    >>> fruitbat.dm_to_redshift(9)
    2.0, 0.0

    >>> fruitbat.dm_to_redshift(9, dm_err=4)
    2.0, 1.0


.. |Travis| image:: https://travis-ci.com/abatten/fruitbat.svg?branch=master
    :target: https://travis-ci.com/abatten/fruitbat

.. |Docs| image:: https://readthedocs.org/projects/fruitbat/badge/?version=latest
    :target: https://fruitbat.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
