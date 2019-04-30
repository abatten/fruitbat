Getting started
===============

Installation
------------
You can install the latest release of **fruitbat** from PyPi_ by running 
the following:

.. code-block:: none

    pip install fruitbat

You can install the latest development version of **fruitbat** by cloning 
the repository_:

.. code-block:: none
    
    git clone https://github.com/abatten/fruitbat
    cd fruitbat
    pip install .

Requirements
************
Below are the listed requirements for running **fruitbat** and the purpose for
each requirement.

 - numpy: Array manipulation

 - scipy: Modules for integration and interpolation

 - astropy: Modules for cosmology, coordinates, constants and units

 - matplotlib: Modules for plotting

 - pandas: Reading csv files from FRBCAT

 - pyymw16: Python wrapper for YMW16 galactic dispersion measure model.

 - e13tools: Utility tools for writing docstrings.


Running Tests
*************
**Fruitbat** uses pytests to perform all its testing. If you would like to run the tests
yourself you will need to install all the packages in `requirements_dev`_. The easiest way of
doing this is by cloning the `source directory`_, install all requirements and running ``pytest``:

.. code-block:: none

    git clone https://github.com/abatten/fruitbat
    cd fruitbat
    pip install -r requirements_dev.txt
    pytest

Pytest should take around 5 minutes.

.. _requirements_dev: https://raw.githubusercontent.com/abatten/fruitbat/master/requirements_dev.txt
.. _source directory: https://github.com/abatten/fruitbat


Simple Example
--------------

A detailed explanation of this example can be viewed at `Using Fruitbat`_.

.. code-block:: python

    >>> import fruitbat

    # Create a Frb Object with DM and Galactic Coordinates
    >>> FRB180110 = fruitbat.Frb(715.7, gl="7.8", gb="-51.9", name="FRB180110")

    # Calculate the DM contribution from the Milky Way
    >>> FRB180110.calc_dm_galaxy()

    # Calculate the Redshift of the FRB using the relation from Zhang (2018)
    >>> FRB180110.calc_redshift(method="Zhang2018", cosmology="Planck18")



.. _repository: https://github.com/abatten/fruitbat
.. _PyPI: https://pypi.org/project/fruitbat
.. _Pyymw16: https://github.com/telegraphic/pyymw16
.. _Using Fruitbat: https://fruitbat.readthedocs.io/en/latest/user_guide/using_fruitbat.html
