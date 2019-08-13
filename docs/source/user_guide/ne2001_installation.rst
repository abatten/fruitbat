NE2001 Installation Instructions
================================

When you install ``fruitbat`` it automatically installs the YMW16 galactic 
electron density model ``pyywm16`` from PyPi. However, if you have the python
port of the NE2001 model by Ben Baror and JXP you can use it with ``fruitbat``
to estimate the galactic contribution to dispersion measure.

The python port that is compatible with ``fruitbat`` is found at: 
https://github.com/FRBs/ne2001. 

You will need to download and install the NE2001 model to use it.


Installing the NE2001 model
---------------------------

.. code-block:: python

    >>> git clone https://github.com/FRBs/ne2001
    >>> cd ne2001
    >>> pip install .

Using the NE2001 model
----------------------

After installing the NE2001 model, it behaves exactly the same as the YNW16 model.
To specify that you want to use the NE2001 model pass `"ne2001"` as a keyword in
:meth:`~fruitbat.Frb.calc_dm_galaxy`.

.. code-block:: python

    >>> import fruitbat
    >>> FRB190523 = fruitbat.Frb(760.8, gl="117.03", gb="44")
    >>> dm_galaxy_ymw16 = FRB190523.calc_dm_galaxy(model="ymw16")
    >>> dm_galaxy_ne2001 = FRB190523.calc_dm_galaxy(model="ne2001")
    >>> print(dm_galaxy_ymw16, dm_galaxy_ne2001)
    29.88017464 pc / cm3    36.87013932 pc / cm3

