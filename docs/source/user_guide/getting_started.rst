Getting started
===============

Installation
------------
**fruitbat** can be easily installed by either cloning the `repository`_ and 
installing it manually::

    $ git clone https://github.com/abatten/fruitbat
    $ cd fruitbat
    $ pip install .

or by installing it directly from `PyPI`_ with::

    $ pip install fruitbat

You can import **fruitbat** as a package using :python:``import fruitbat``.


Examples
--------

.. code-block:: python

    import fruitbat

    # Create Frb Object
    FRB121110 = fruitbat.Frb("FRB121110", dm=534, dm_galaxy=30)

    # Calculate the redshift using Inoue 2004 with Planck 2018 cosmology
    FRB121110.calc_redshift(method='inoue2004', cosmology='planck2018+bao')


.. _repository: https://github.com/abatten/fruitbat
.. _PyPI: https://pypi.org/project/fruitbat
