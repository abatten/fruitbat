"""
Module for defining different cosmologies
"""
from __future__ import print_function, absolute_import, division

from astropy import units as u
import astropy.cosmology
from astropy.cosmology.core import (FlatLambdaCDM, FlatwCDM, LambdaCDM, wCDM)
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


def create_cosmology(parameters=None, name=None):
    """
    A wrapper to create custom astropy cosmologies.

    The only avaliable cosmology types in this method are: FlatLambdaCDM,
    FlatwCDM, LambdaCDM and wCDM. See `astropy.cosmology`_ for more details on
    these types of cosmologies. To create a cosmology of a type that isn't
    listed above, it will have to be created directly using astropy.cosmology.

    Parameters
    ----------
    parameters: dict or None
        A dictionary containing the cosmological parameters. The names of the
        parameters must conform to the same format as the parameters used in
        astropy.cosmology. If `parameters` is *None* then default values for
        each parameter is used.

    name: str or None, optional
        The name of the cosmology. Default: *None*

    Returns
    -------
    cosmology


    Notes
    -----
    Default parameter values:

    .. code-block:: python

        params = {'H0': 70, 'Om0': 0.3, 'Oc0': 0.26, 'Ob0': 0.04, 'Neff': 3.04,
                  'flat': True, 'Tcmb0': 0.0, 'm_nu': 0.0, 'w0': -1}

    If ``'flat'`` is set to ``False`` then a value of ``'Ode0'`` (current dark
    energy density) should be specified.

    .. _astropy.cosmology: http://docs.astropy.org/en/latest/cosmology/index.html

    """

    # First check if the cosmology is flat or curved. Then check to see if dark
    # energy is parameterised by a cosmolological constant or an equation of
    # state.

    # Set the default parameters:
    params = {'H0': 70, 'Om0': 0.3, 'Oc0': 0.26, 'Ob0': 0.04, 'w0': -1,
              'Neff': 3.04, 'flat': True, 'Tcmb0': 0.0, 'm_nu': 0.0}

    # Override default parameters with supplied parameters
    if parameters is not None:
        params.update(parameters)

    if params["flat"]:
        if params['w0'] is not -1:
            cosmo = FlatwCDM(H0=params['H0'], Om0=params['Om0'],
                             w0=params['w0'], Tcmb0=params['Tcmb0'],
                             Neff=params['Neff'], Ob0=params['Ob0'],
                             m_nu=u.Quantity(params['m_nu'], u.eV), name=name)

        else:
            cosmo = FlatLambdaCDM(H0=params['H0'], Om0=params['Om0'],
                                  Tcmb0=params['Tcmb0'], Neff=params['Neff'],
                                  Ob0=params['Ob0'], name=name,
                                  m_nu=u.Quantity(params['m_nu'], u.eV))

    else:
        if params['w0'] is not -1:
            cosmo = wCDM(H0=params['H0'], Om0=params['Om0'],
                         Ode0=params['Ode0'], w0=params['w0'],
                         Tcmb0=params['Tcmb0'], Neff=params['Neff'],
                         m_nu=u.Quantity(params['m_nu'], u.eV), name=name,
                         Ob0=params['Ob0'])

        else:
            cosmo = LambdaCDM(H0=params['H0'], Om0=params['Om0'],
                              Ode0=params['Ode0'], Tcmb0=params['Tcmb0'],
                              Neff=params['Neff'], Ob0=params['Ob0'],
                              m_nu=u.Quantity(params['m_nu'], u.eV), name=name)
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
    reference=("Planck 2018 results. VI. Cosmological Parameters, A&A,"
               "submitted, Table 2 (TT, TE, EE + lowE + lensing + BAO)")
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
        "Planck18": Planck18(),
        "EAGLE": Planck13(),
    }
    return cosmologies


def keys():
    """
    Returns a string constaining all the keywords for builtin cosmologies.
    """
    cosmologies = builtin()
    keys = ", ".join(cosmologies.keys())
    return keys
