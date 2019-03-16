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
:attr:`~fruitbat.Frb.name`, which can be useful when calculating redhsifts of
many FRBs at once.

.. code-block:: python
    
    >>> import fruitbat
    >>> frb = fruitbat.Frb(635.1, name="FRB190229")
   

Redshift Estimation
*******************

.. code-block:: python
    
    >>> frb = fruitbat.Frb(635.1)
    >>> frb.calc_redshift()
    <Quantity 0.63199287>

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
a value of :attr:`~fruitbat.Frb.dm_galaxy`. This 

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

Cosmology
~~~~~~~~~
