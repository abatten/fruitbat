---
title: 'Fruitbat: A tool for estimating distances to FRBs'
tags:
  - astrophysics
  - frb
  - radio
  - transients
authors:
  - name: Adam J. Batten
    orcid: 0000-0001-7559-6488
    affiliation: "1, 2"
affiliations:
- name: Centre for Astrophysics and Supercomputing, Swinburne University of Technology, PO Box 218, Hawthorn, VIC 3122, Australia
  index: 1
- name: ARC Centre of Excellence for All Sky Astrophysics in 3 Dimensions (ASTRO 3D)
  index: 2
date: 1 April 2019
bibliography: paper.bib
---

# Summary

``fruitbat`` is a Python 3 package for estimating redshifts, energies and the 
galactic dispersion measure (DM) contributions of fast radio bursts (FRBs).

FRBs are a class of short duration (~ 30 ms) transient radio sources with an 
unknown extragalactic origin (@Lorimer2007, @Thornton2013, @Petroff2015, 
@CHIME2019). There is currently at least 65 confirmed FRB detections 
(@Petroff2016), with the new CHIME telescope expected to detect of the order 
10 FRBs per day (@Chawla2017).

The defining feature of FRBs that sets them apart from other radio transient 
events is their integrated free electron column density, called DM. FRBs have 
DM values significantly higher than the estimated Milky Way contribution along
the line-of-sight, leading many to believe that their origin is extragalactic.
The extragalactic origin of FRBs was confirmed with the host galaxy 
localisation of a repeating FRB (FRB 121102) to a redshift of $z = 0.197$ 
(@Tendulkar2017).

However, most telescopes do not have the capability to localise host galaxies 
of FRBs and the upper limit of their redshifts must be estimated using a 
DM-redshift relation; which is typically calculated via analytical means 
(e.g. @Ioka2003, @Inoue2004, @Zhang2018).

We have built ``fruitbat`` as a tool in assisting the estimating of redshifts 
and galactic DM values of FRBs. ``fruitbat`` generates and utilises 'look-up 
tables' (saved instance methods of scipy.interp1d) of existing DM-redshift 
relations found in the literature (@Ioka2003, @Inoue2004, @Zhang2018) in 
conjunction with cosmological parameters determined from both the WMAP and 
Planck missions. 

This feature of ``fruitbat`` allows the user to independantly choose the 
DM-redshift relation and the cosmological parameters which was typically not 
an option when using the relations from the literature. Additionally, 
``fruitbat`` explicitly integrates the entire DM-z relation at each redshift 
instead of assuming an average value for across redshifts which introduces a 6%
error (see Equation (6) in @Zhang2018). ``fruitbat`` has been built to allow 
for future non-analytical DM-redshift relations such as those from cosmological
    simulations (Batten et al. in prep.).

To account for the galactic DM contribution due to electrons the the 
interstellar and circumgalatic medium, ``fruitbat`` utilises the YMW16 galactic
free electron density model (@Yao2017) to estimate the line-of-sight DM of the
Milky Way.

``fruitbat`` additionally has the feature to estimate the following FRB values:
* Burst energy
* Average luminosity 
* Comoving distance
* Luminosity distance
* Fluence

``fruitbat`` has been used in Price et al. submitted to estimate the redshift 
of FRB 180301.

``fruitbat`` is released under the BSD 3-Clause licence, and is avaliable from
PyPi via ``pip``; source code for ``fruitbat`` can be found at 
https://github.com/abatten/fruitbat.


# Acknowledgements
Adam Batten would like to thank Alan Duffy (@astroduff), Ellert van der Velden 
(@1313e) and Daniel Price (@telegraphic).

# References
