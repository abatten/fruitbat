|Travis|

FRBz
====

FRBz is a package estimate the redshift of a Fast Radio Burst (FRB) from its
dispersion measure. 

Installation
------------

Run the following to install::

    pip install frbz

Usage
-----

The main function you want to use is `dm_to_redshift()`::

    >>> import frbz
    >>> frbz.dm_to_redshift(9)
    2.0, 0.0

    >>> frbs.dm_to_redshift(9, dm_err=4)
    2.0, 1.0


.. |Travis| image:: https://travis-ci.com/abatten/frbz.svg?token=cSfgUVgVHZsxUNLefqMs&branch=master
    :target: https://travis-ci.com/abatten/frbz
