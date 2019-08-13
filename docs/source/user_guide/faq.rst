FAQ
===

Why is it called fruitbat?
--------------------------
A lot of time was spent trying to figure out a 'nice' acronym, however none of the
potential names seemed to fit. So I decided to just name it after an animal. 
After searching through lists of Australian animals I remembered the family of
fruitbats (grey-headed flying foxes) that live in the tree in the front of my house.
Fruitbat had the letters FRB in order and their orange necks reminded me of waterfall
plots so I decided to go with the name Fruitbat.

I'd rather use the NE2001 model. Is this compatible with Fruitbat?
------------------------------------------------------------------
Yes! However you will need to install the python port of the NE2001 model.
You can find installation and usage instructions here_.

Alternatively, you can calculate the galactic DM using any version of the NE2001 model
first and then use that value for :attr:`~fruitbat.Frb.dm_galaxy`.

.. code-block:: python

    >>> ne2001_dm = 34.5
    >>> import fruitbat
    >>> FRB = fruitbat.Frb(dm=325, dm_galaxy=ne2001_dm, gl="32.4", gb="16.2")


.. _here: https://fruitbat.io/en/latest/user_guide/ne2001_installation.html
