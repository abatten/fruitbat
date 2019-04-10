Using Fruitbat
==============

Here is the guide of how to use the **fruitbat** package.

Frb Class
~~~~~~~~~
Most calculations in **fruitbat** are centred around the class 
:class:`~fruitbat.Frb`. Hence you will need to define a :class:`fruitbat.Frb` 
object for each FRB. The minimum required to define a :class:`~fruitbat.Frb` 
is the observed dispersion measure (DM). 

.. code-block:: python

    >>> import fruitbat
    >>> frb = fruitbat.Frb(635.1)

It's also possible to 'name' the FRB using the keyword 
:attr:`~fruitbat.Frb.name`, which can be useful when calculating redshifts of
many FRBs at once.

.. code-block:: python
    
    >>> import fruitbat
    >>> frb = fruitbat.Frb(635.1, name="FRB190229")
   

Redshift Estimation
*******************
To estimate the redshift of the FRB, use the method 
:meth:`~fruitbat.Frb.calc_redshift`. Unless otherwise specified this assumes
that the entire DM contribution is due to the IGM. To see how to account for the 
Milky Way and host galaxy contributions to the DM see the sections on 
`Galactic Dispersion Measure`_ and `Host Galaxy Dispersion Measure`_
respectively.

.. code-block:: python
    
    >>> frb = fruitbat.Frb(635.1)
    >>> frb.calc_redshift()
    <Quantity 0.63199287>

The :meth:`~fruitbat.Frb.calc_redshift` also has keywords to select alternative
IGM models and cosmologies when estimating the redshift of the FRB (Default
method and cosmology are ``Inoue2004`` and ``Planck2018`` respectively).

.. code-block:: python

    >>> frb.calc_redshift(method="Zhang2018")
    <Quantity 0.70986024>

    >>> frb.calc_redshift(method="Ioka2003", cosmology="Planck13")
    <Quantity 0.52776778>


Currently avaliable methods in **fruitbat** include: ``Ioka2003``, 
``Inoue2004``, ``Zhang2018``. Currently avaliable cosmologies include: 
``WMAP5``, ``WMAP7``, ``WMAP9``, ``Planck13``, ``Planck15``, ``Planck18``, 
``EAGLE``. It should be of note that ``EAGLE`` uses the ``Planck13`` cosmology
but is listed here for convenience.

.. _Galactic Dispersion Measure:

Galactic Dispersion Measure
...........................

When estimating the redshift of an FRB, it is necessary to account for the 
Milky Way's contribution to the observed DM. Depending on the location on the 
sky, this contribution can be as low as :math:`~30\ \rm{pc\ cm^{-3}}` out of 
the disk of the Milky Way and exceeding :math:`1000\ \rm{pc\ cm^{-3}}` through
the disk. Without accuretly accounting for this contribution, redshift
estimates of FRB would be significantly higher than their physical redshift.

Within **fruitbat** there are three
main ways to account for the galactic contribution: 
:meth:`~fruitbat.Frb.calc_dm_galaxy`, :attr:`~fruitbat.Frb.dm_galaxy` or
:attr:`~fruitbat.Frb.dm_excess`.  

Method 1: calc_dm_galaxy()
--------------------------
The first and easiest way to account for the galactic contribution is to
provide the sky coordinates of the FRB when instantiating the object, then 
call :meth:`~fruitbat.Frb.calc_dm_galaxy`. The 
:meth:`~fruitbat.Frb.calc_dm_galaxy` method of :class:`~fruitbat.Frb` estimates
the total DM contribution due to the Milky Way along the line of sight of the 
FRB using the YMW16 galactic free electron model. 



.. code-block:: python

    >>> frb = fruitbat.Frb(635.1, gl="35.1", gb="12.5")
    >>> frb.calc_dm_galaxy()
    <Quantity 114.27922821 pc / cm3>

    >>> frb = fruitbat.Frb(635.1, raj="18:10:34.8668", decj="7:33:35.9289")
    >>> frb.calc_dm_galaxy()
    <Quantity 114.27922821 pc / cm3>
    

The sky coordinates can be in either ICRS or Galactic units. The 
:meth:`~fruitbat.Frb.calc_dm_galaxy` method will calculate the 
:attr:`~fruitbat.Frb.dm_excess` by subtracting the estimated 
:attr:`~fruitbat.Frb.dm_galaxy` from the observed DM. After calculating 
:attr:`~fruitbat.Frb.dm_galaxy`, calling :meth:`~fruitbat.Frb.calc_redshift`
will automatically use the calculated :attr:`~fruitbat.Frb.dm_excess` to 
estimate the redshift.

.. code-block:: python

    >>> frb.calc_redshift()
    <Quantity 0.52244866>    

Method 2: dm_galaxy
-------------------
The second method to account for the galactic dispersion meausre is to provide
a value of :attr:`~fruitbat.Frb.dm_galaxy`. 

.. code-block:: python

    >>> frb = fruitbat.Frb(635.1, dm_galaxy=114.28)
    >>> frb.calc_redshift()
    <Quantity 0.52244791>


Method 3: dm_excess
-------------------
The third and final method is to directly specify the 
:attr:`~fruitbat.Frb.dm_excess`.

.. code-block:: python

    >>> frb = fruitbat.Frb(635.1, dm_excess=520.82)
    >>> frb.calc_redshift()
    <Quantity 0.52244791>


.. _Host Galaxy Dispersion Measure:

Host Galaxy Dispersion Measure
..............................

Calculating Distances
*********************

Luminosity Distance
...................

.. code-block:: python

    >>> frb = fruitbat.Frb(635.1, gl="35.1", gb="12.5")
    >>> frb.calc_dm_galaxy()
    <Quantity 114.27922821 pc / cm3>

    >>> frb.calc_redshift()
    <Quantity 0.52244866>

    >>> frb.calc_luminosity_distance()
    <Quantity 3075.79950018 Mpc>

Comoving Distance
.................

.. code-block:: python

    >>> frb.calc_comoving_distance()
    <Quantity 2020.29768846 Mpc>


Calculating Energy
******************

.. code-block:: python

    >>> frb.fluence = 2.0
    >>> frb.obs_freq_central = 1600
    >>> frb.calc_energy()
    <Quantity 2.37921847e+40 erg>


Custom Lookup Tables
~~~~~~~~~~~~~~~~~~~~

Custom Methods
**************

.. code-block:: python

    >>> def simple_dm(z):
        dm = 1200 * z
        return dm
    >>> fruitbat.add_method("simple_dm", simple_dm)
    ['Ioka2003', 'Inoue2004', 'Zhang2018', 'simple_dm']



Custom Cosmologies
******************

.. code-block:: python

    >>> params = {"H0": 72.4, "Om0": 0.26}
    >>> new_cosmology = fruitbat.cosmology.create_cosmology(parameters=params)
    >>> fruitbat.add_cosmology("new_cosmology", new_cosmology)
    >>> fruitbat.avaliable_cosmologies()
    ['WMAP5', 'WMAP7', 'WMAP9', 'Planck13', 'Planck15', 'Planck18', 'EAGLE', 'new_cosmology']


Custom Table
************

.. code-block:: python

    >>> def simple_dm(z):
        dm = 1200 * z
        return dm
    >>> fruitbat.add_method("simple_dm", simple_dm)
    >>> fruitbat.table.create("simple_dm")
    >>> frb = fruitbat.Frb(1200)
    >>> frb.calc_redshift(method="simple_dm")
    <Quantity 1.>

