"""
Module for defining different cosmologies
"""
from astropy import units as u


def avaliable_cosmologies():
    """
    Returns
    -------
    dict
        A dictionary containing the parameters for each cosmology.
    """
    cosmology_dict = {
        "wmap2013": wmap2013(),
        "planck2013": planck2013(),
        "planck2015": planck2015(),
        "planck2018": planck2018(),
        "eagle": eagle()
    }

    return cosmology_dict 


def cosmology_keys():
    """
    Returns
    -------
    list
        A list containing the keyword for all avaliable cosmologies.
    """

    return avaliable_cosmologies().keys()


def planck2018():
    """

    https://ui.adsabs.harvard.edu/#abs/arXiv:1807.06209

    Returns
    -------
    dict
        Planck 2018 best fit with cmb lensing and bao
    """
    
    planck2018 = {}

    planck2018["HO"] = 67.66 * u.km / u.s / u.Mpc
    planck2018["Omega_m"] = 0.3111
    planck2018["Omega_L"] = 0.6889
    planck2018["Omega_b"] = 0.04897
    planck2018["Omega_c"] = 0.2606
    planck2018["sigma8"] = 0.8102 * u.Mpc

    return planck2018


def planck2015():
    """
    Returns
    -------
    dict
        Planck 2015 best fit + CMB lensing + baryonic acoustic oscillations
    """
    
    planck2015 = {}

    planck2015["HO"] = 67.74 * u.km / u.s / u.Mpc 
    planck2015["Omega_m"] = 0.3089
    planck2015["Omega_L"] = 0.6911
    planck2015["Omega_b"] = 0.04860
    planck2015["Omega_c"] = 0.2589
    planck2015["sigma8"] = 0.8159 * u.Mpc

    return planck2015


def planck2013():
    """
    Returns
    -------
    dict
        Planck 2013 best fit + bao
    """
    
    planck2013 = {}

    planck2013["HO"] = 67.80 * u.km / u.s / u.Mpc 
    planck2013["Omega_m"] = 0.3063
    planck2013["Omega_L"] = 0.692
    planck2013["Omega_b"] = 0.04816
    planck2013["Omega_c"] = 0.2582
    planck2013["sigma8"] = 0.826 * u.Mpc

    return planck2013


def wmap2013():
    """
    Returns
    -------
    dict
        wmap 2013 best fit final data release
    """
    
    wmap2013 = {}

    wmap2013["HO"] = 69.32 * u.km / u.s / u.Mpc 
    wmap2013["Omega_m"] = 0.2865
    wmap2013["Omega_L"] = 0.7135
    wmap2013["Omega_b"] = 0.04628
    wmap2013["Omega_c"] = 0.2402
    wmap2013["sigma8"] = 0.820 * u.Mpc

    return wmap2013

def eagle():
    """
    Returns
    -------
    dict
        The cosmology that is used in the EAGLE simulations
    """

    eagle = {}

    eagle["HO"] = 67.77 * u.km / u.s / u.Mpc
    eagle["Omega_m"] = 0.307
    eagle["Omega_L"] = 0.693
    eagle["Omega_b"] = 0.04825
    eagle["sigma8"] = 0.8288
    eagle["ns"] = 0.9611
    eagle["Y"] = 0.248

    return eagle
