Getting started
===============

Installation
------------
You can install the latest release of ``fruitbat`` from PyPi_ by running 
the following::

    pip install fruitbat

You can install the latest deveopment version of ``fruitbat`` by cloning 
the repository::
    
    git clone https://github.com/abatten/fruitbat
    cd fruitbat
    pip install .

Requirements
************
* Numpy >= 1.12.0
* Astropy >= 2.0.0
* Scipy >= 1.0.0
* Pyymw16 >= 2.0.4

Pyymw16_ is a python wrapper for the YMW16 galactic dispersion measure model.


Using Fruitbat
--------------

A detailed explanation of this example can be viewed at `Using Fruitbat`_.


Example Calculation
*******************

.. code-block:: python

    import fruitbat

    # Create a Frb Object with DM and Galactic Coordinates
    FRB180110 = fruitbat.Frb("FRB180110", dm=715.7, gl="7.8", gb="-51.9")

    # Calculate the DM contribution from the Milky Way
    FRB180110.calc_dm_galaxy(model="ymw16")

    # Calculate the Redshift of the FRB using the relation from Zhang (2018)
    FRB180110.calc_redshift(method="zhang2018", cosmology="planck2018")



.. _repository: https://github.com/abatten/fruitbat
.. _PyPI: https://pypi.org/project/fruitbat
.. _Pyymw16: https://github.com/telegraphic/pyymw16
.. _Using Fruitbat: https://fruitbat.readthedocs.io/en/latest/user_guide/using_fruitbat.html
