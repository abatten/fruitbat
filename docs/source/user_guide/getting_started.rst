Getting started
===============

Installation
------------
You can install the latest release of **fruitbat** from PyPi_ by running 
the following:

.. code-block:: none

    pip install fruitbat

.. _PyPi: https://pypi.python.org/pypi/fruitbat 

You can install the latest development version of **fruitbat** by cloning 
the repository_:

.. code-block:: none
    
    git clone https://github.com/abatten/fruitbat
    cd fruitbat
    pip install .

.. _repository: https://github.com/abatten/fruitbat

If you are installing the latest development version of ``fruitbat`` then you 
will also need to install git-lfs. Instructions for installing git-lfs for
your operating system can be found here_.

.. _here: https://help.github.com/en/articles/installing-git-large-file-storage

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


Simple Example
--------------

A detailed explanation of this example can be viewed at `Using Fruitbat`_.

.. code-block:: python

    import fruitbat

    # Create a Frb Object with DM and Galactic Coordinates
    FRB180110 = fruitbat.Frb(715.7, gl="7.8", gb="-51.9", name="FRB180110")

    # Calculate the DM contribution from the Milky Way
    FRB180110.calc_dm_galaxy()

    # Calculate the Redshift of the FRB using the relation from Zhang (2018)
    FRB180110.calc_redshift(method="Zhang2018", cosmology="Planck18")



.. _Using Fruitbat: https://fruitbat.readthedocs.io/en/latest/user_guide/using_fruitbat.html
