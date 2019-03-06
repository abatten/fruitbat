"""
Module for defining different cosmologies
"""
from astropy import units as u
import astropy.cosmology
from astropy.cosmology.core import FlatLambdaCDM
from astropy.cosmology.core import FlatwCDM
from astropy.cosmology.core import LambdaCDM
from astropy.cosmology.core import wCDM
from e13tools import docstring_copy

__all__ = ["WMAP5", "WMAP7", "WMAP9", "Planck13", "Planck15", "Planck18",
          "create_cosmology", "builtin", "keys"]

@docstring_copy(astropy.cosmology.WMAP9)
def WMAP9():
    return astropy.cosmology.WMAP9

@docstring_copy(astropy.cosmology.WMAP7)
def WMAP7():
    return astropy.cosmology.WMAP7

@docstring_copy(astropy.cosmology.WMAP5)
def WMAP5():
    return astropy.cosmology.WMAP5

@docstring_copy(astropy.cosmology.Planck13)
def Planck13():
    return astropy.cosmology.Planck13

@docstring_copy(astropy.cosmology.Planck15)
def Planck15():
    return astropy.cosmology.Planck15

def Planck18():
    """
    Planck18 instance of FlatLambdaCDM cosmology

    (from Planck 2018 results. VI. Cosmological Parameters, 
    A&A, submitted, Table 2 (TT, TE, EE + lowE + lensing + BAO))
    """
    cosmo = create_cosmology(name="Planck18", parameters=planck18)    
    return cosmo


def create_cosmology(name, parameters=None, default=False):
    """
    A wrapper to create custom astropy cosmologies. 
    
    The only avaliable cosmology types in this method are: FlatLambdaCDM, 
    FlatwCDM, LambdaCDM and wCDM. See `astropy.cosmology`_ for more details on 
    these types of cosmologies. To create a cosmology of a type that isn't 
    listed above, it will have to be created directly using astropy.cosmology.

    Parameters
    ----------
    name: str
        The name of the cosmology. 

    parameters: dict or None
        A dictionary containing the cosmological parameters. The names of the 
        parameters must conform to the same format as those in 
        astropy.cosmology. The minimum parameters needed to specify
        a cosmology is: 'H0', 'Om0', 'Tcmb0', 'Neff', 'm_nu', 'Ob0'

    default: bool, optional
        Create a cosmology with default values for parameters. Should only be
        used for testing. Default: False

    .. _astropy.cosmology: http://docs.astropy.org/en/latest/cosmology/index.html

    """

    # First check if the cosmology is flat or curved. Then check to see if dark
    # energy is parameterised by a cosmolological constant or an equation of 
    # state.

    if parameters is not None:
        if parameters["flat"]:
            if 'w0' in parameters:
                cosmo = FlatwCDM(H0=parameters['H0'],
                                 Om0=parameters['Om0'],
                                 w0=parameters['w0'],
                                 Tcmb0=parameters['Tcmb0'],
                                 Neff=parameters['Neff'],
                                 m_nu=u.Quantity(parameters['m_nu'], u.eV),
                                 name=name, 
                                 Ob0=parameters['Ob0'])
                cosmotype = "FlatwCDM"

            else:
                cosmo = FlatLambdaCDM(H0=parameters['H0'],
                                      Om0=parameters['Om0'],
                                      Tcmb0=parameters['Tcmb0'],
                                      Neff=parameters['Neff'],
                                      m_nu=u.Quantity(parameters['m_nu'], u.eV),
                                      name=name, 
                                      Ob0=parameters['Ob0'])
                cosmotype = "FlatLambdaCDM"

        else:
            if 'w0' in parameters:
                cosmo = wCDM(H0=parameters['H0'],
                                  Om0=parameters['Om0'],
                                  Ode0=parameters['Ode0'],
                                  w0=parameters['w0'],
                                  Tcmb0=parameters['Tcmb0'],
                                  Neff=parameters['Neff'],
                                  m_nu=u.Quantity(parameters['m_nu'], u.eV),
                                  name=name, 
                                  Ob0=parameters['Ob0'])
                cosmotype = "wCDM"

            else:
                cosmo = LambdaCDM(H0=parameters['H0'],
                                  Om0=parameters['Om0'],
                                  Ode0=parameters['Ode0'],
                                  Tcmb0=parameters['Tcmb0'],
                                  Neff=parameters['Neff'],
                                  m_nu=u.Quantity(parameters['m_nu'], u.eV),
                                  name=name, 
                                  Ob0=parameters['Ob0'])
                cosmotype = "LambdaCDM"

        if parameters["reference"]:
            docstr =  "{} instance of {} cosmology \n\n(from {}) "
            cosmo.__doc__ = docstr.format(name, cosmotype, 
                                          parameters["reference"])

        else:
            docstr =  "{} instance of {} cosmology"
            cosmo.__doc__ = docstr.format(name, cosmotype)


    # Have the option to create a mock FlatLambdaCDM without a param dict    
    elif default:
            cosmo = FlatLambdaCDM(H0=70,
                                  Om0=0.3,
                                  Tcmb0=0,
                                  Neff=3.04,
                                  m_nu=u.Quantity(0.0, u.eV),
                                  name=name, 
                                  Ob0=0.04)
    else:
        raise ValueError("""Either a parameter dictionary must be supplied"""
                         """ or specify that default=True""")

    return cosmo
            

# Planck 2018 paper VI Table 2 Final column (68% confidence interval)
planck18 = dict(
    Oc0=0.2607,
    Ob0=0.04897,
    Om0=0.3111,
    H0=67.66,
    n=0.9665,
    sigma8=0.8102,
    tau=0.0561,
    z_reion=7.82,
    t0=13.787,
    Tcmb0=2.7255,
    Neff=3.046,
    flat=True,
    m_nu=[0., 0., 0.06],
    reference=("Planck 2018 results. VI. Cosmological Parameters, A&A, submitted,"
               " Table 2 (TT, TE, EE + lowE + lensing + BAO)")
)


def builtin():
    """
    Create a dictionary of the builtin cosmologies with keywords and functions.

    Returns
    -------
    dict
        A dictionary containing the keywords and function for each cosmology.
    """
    cosmologies = {
        "WMAP5": WMAP5(),
        "WMAP7": WMAP7(),
        "WMAP9": WMAP9(),
        "Planck13": Planck13(),
        "Planck15": Planck15(),
        "Planck18": Planck18()
    }

    return cosmologies


def keys():
    """
    Returns the keywords for builtin cosmologies.
    """
    cosmologies = builtin()
    keys = ", ".join(cosmologies.keys())
    return keys


#def planck2018():
#    """
#
#    https://ui.adsabs.harvard.edu/#abs/arXiv:1807.06209
#
#    Returns
#    -------
#    dict
#        Planck 2018 best fit with cmb lensing and bao
#    """
#    
#    planck2018 = {}
#
#    planck2018["HO"] = 67.66 * u.km / u.s / u.Mpc
#    planck2018["Omega_m"] = 0.3111
#    planck2018["Omega_L"] = 0.6889
#    planck2018["Omega_b"] = 0.04897
#    planck2018["Omega_c"] = 0.2606
#    planck2018["sigma8"] = 0.8102 * u.Mpc
#
#    return planck2018
#
#
#def planck2015():
#    """
#    Returns
#    -------
#    dict
#        Planck 2015 best fit + CMB lensing + baryonic acoustic oscillations
#    """
#    
#    planck2015 = {}
#
#    planck2015["HO"] = 67.74 * u.km / u.s / u.Mpc 
#    planck2015["Omega_m"] = 0.3089
#    planck2015["Omega_L"] = 0.6911
#    planck2015["Omega_b"] = 0.04860
#    planck2015["Omega_c"] = 0.2589
#    planck2015["sigma8"] = 0.8159 * u.Mpc
#
#    return planck2015
#
#
#def planck2013():
#    """
#    Returns
#    -------
#    dict
#        Planck 2013 best fit + bao
#    """
#    
#    planck2013 = {}
#
#    planck2013["HO"] = 67.80 * u.km / u.s / u.Mpc 
#    planck2013["Omega_m"] = 0.3063
#    planck2013["Omega_L"] = 0.692
#    planck2013["Omega_b"] = 0.04816
#    planck2013["Omega_c"] = 0.2582
#    planck2013["sigma8"] = 0.826 * u.Mpc
#
#    return planck2013
#
#
#def wmap2013():
#    """
#    Returns
#    -------
#    dict
#        wmap 2013 best fit final data release
#    """
#    
#    wmap2013 = {}
#
#    wmap2013["HO"] = 69.32 * u.km / u.s / u.Mpc 
#    wmap2013["Omega_m"] = 0.2865
#    wmap2013["Omega_L"] = 0.7135
#    wmap2013["Omega_b"] = 0.04628
#    wmap2013["Omega_c"] = 0.2402
#    wmap2013["sigma8"] = 0.820 * u.Mpc
#
#    return wmap2013
#
#def eagle():
#    """
#    Returns
#    -------
#    dict
#        The cosmology that is used in the EAGLE simulations
#    """
#
#    eagle = {}
#
#    eagle["HO"] = 67.77 * u.km / u.s / u.Mpc
#    eagle["Omega_m"] = 0.307
#    eagle["Omega_L"] = 0.693
#    eagle["Omega_b"] = 0.04825
#    eagle["sigma8"] = 0.8288
#    eagle["ns"] = 0.9611
#    eagle["Y"] = 0.248
#
#    return eagle
