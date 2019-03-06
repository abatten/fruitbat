from e13tools import docstring_substitute

from . import utils
from ._fruitbatstrings import dm_units_doc
from .cosmology import keys as cosmo_keys, builtin

__all__ = ["redshift", "methods"]


def methods(string=False):
    """
    Defines the list of avaliable method keywords.

    Methods currently avaliable: ioka2003, inoue2004, zhang2018

    Parameters
    ----------
    string: bool, optional
        If True, return a string of keywords instead of a list.

    Returns
    -------
    list or str:
        A list containing the valid method keywords. If ``string=True`` it
        returns a single string listing all the keywords.
    """
    methods = ["ioka2003", "inoue2004", "zhang2018"]

    if string:
        methods = ", ".join(methods)

    return methods


@docstring_substitute(dmunits=dm_units_doc, methods=methods(string=True), 
                      cosmo=cosmo_keys())
def redshift(dm, dm_uncert=0.0, method='inoue2004', cosmology='Planck18'):
    """
    Returns the redshift of a given dispersion measure using a
    specified DM-z relation.

    Parameters
    ----------
    dm : float
        Dispersion Measure. Units: %(dmunits)s

    dm_uncert : float or None
        The uncertainty in the dispersion measure. Units: %(dmunits)s

    method : string, optional
        The DM-z relation to use to calculate the redshift.
        Avaliable methods are: %(methods)s. Default: `'inoue2004'`

    cosmology : string, optional
        Avaliable cosmologies: %(cosmo)s. Default: `'planck2018'`

    Returns
    -------
    z : float
        Redshift

    z_err : float
        The uncertianty in the redshift estimation. If `dm_uncert` is `None`
        then `z_err` = 0.

    Notes
    -----

    Cosmology_ has a list of the cosmological parameters used in each
    cosmology method.

    .. _Cosmology: https://fruitbat.readthedocs.io/en/latest/cosmology.html
    """

    valid_methods = methods()

    if method not in valid_methods:
        raise ValueError("""Method '{}' is not a valid method.
            Valid methods are: {}""".format(method, methods(string=True)))

    if cosmology not in builtin().keys():
        raise ValueError("""Cosmology '{}' is not a valid cosmology.
            Valid cosmologies are: {}""".format(cosmology, cosmo_keys()))

    z = _get_redshift_from_table(dm, method, cosmology)

    return z


#def _redshift_batten2019(dm, dm_uncert=0.0):
#    """
#    Calculates a redshft from a dispersion measure using the DM-z
#    relation from Batten, A. J. 2019, ....
#    """
#    return 12.0
#
#@_docstr_sub(dm_units=_dm_units_doc)
#def _redshift_inoue2004(dm, dm_uncert=0.0, cosmology="planck2018"):
#    """
#    Calculates a redshift from a dispersion measure using the DM-z
#    relation from Inoue2004_
#
#    Parameters
#    ----------
#    dm : float
#        Units: %(dm_units)s
#
#    dm_uncert : float, optional
#        Default: :0.0
#
#    cosmology : str, optional
#        Default: `'planck2018'`
#
#    Returns
#    -------
#    z : float
#        The redshift of the FRB.
#    z_uncert : tuple of floats
#        The uncertainty in the redshift of the FRB.
#
#    Notes
#    -----
#    Inoue2004_ presents the following dispersion measure and redshift relation.
#   
#    .. math ::
#        \\mathrm{DM}(z) =0.92\\times 10^{-5} \\Omega_b h^2 c \\int_0^z
#        \\frac{1 + z'}{H_0\\left(\\Omega_m (1 + z)^3 +
#        \\Omega_\\Lambda\\right)^{1/2}} dz'
#
#
#    .. _Inoue2004: http://adsabs.harvard.edu/abs/2004MNRAS.348..999I
#    """
#
#    cosmo_dict = {
#        "wmap2013": "inoue2004_wmap2013.npy",
#        "planck2013": "inoue2004_planck2013.npy",
#        "planck2015": "inoue2004_planck2015.npy", 
#        "planck2018": "inoue2004_planck2018.npy",
#        "eagle": "inoue2004_eagle.npy",
#    }
#
#    lookup_table = utils.load_lookup_table(cosmo_dict[cosmology])
#
#    z = _get_redshift_from_table(lookup_table, dm, dm_uncert) 
#    #z, dz_low, dz_high = _get_redshift_from_table(lookup_table, dm, dm_uncert) 
#    return z#, (dz_low, dz_high)
#
#
#def _redshift_ioka2003(dm, dm_uncert=0.0):
#    """
#    Calculates a redshift from a dispersion measure using the DM-z
#    relation from Ioka, K. 2003, ApJL, 598, L79
#    """
#    cosmo_dict = {
#        "wmap2013": "ioka2003_wmap2013.npy",
#        "planck2013": "ioka2003_planck2013.npy",
#        "planck2015": "ioka2003_planck2015.npy",
#        "planck2018": "ioka2003_planck2018.npy",
#        "eagle": "ioka2003_eagle.npy"
#        }
#    lookup_table = utils.load_lookup_table(cosmo_dict[cosmology])
#    z = _get_redshift_from_table(lookup_table, dm, dm_uncert)
#
#    return z
#
#
#def _redshift_zhang2018(dm, dm_uncert=0.0, cosmology="planck2018"):
#    """
#    Calculates a redshift from a dispersion measure using the DM-z
#    relation from Zhang2018_
#
#
#    .. _Zhang2018: https://ui.adsabs.harvard.edu/#abs/arXiv:1808.05277
#    """
#
#    cosmo_dict = {
#        "wmap2013": "zhang2018_wmap2013.npy",
#        "planck2013": "zhang2018_planck2013.npy",
#        "planck2015": "zhang2018_planck2015.npy",
#        "planck2018": "zhang2018_planck2018.npy",
#        "eagle": "zhang2018_eagle.npy"
#        }
#    lookup_table = utils.load_lookup_table(cosmo_dict[cosmology])
#    z = _get_redshift_from_table(lookup_table, dm, dm_uncert)
#
#    return z
#

def _get_redshift_from_table(dm, method, cosmology):
    """
    Loads a lookup table for a corresponding method and cosmology and reads
    the value of the redshift correspoding the the DM value.

    Parameters
    ----------
    dm : int or float

    method: str

    cosmology: str
        The cosmology keyword 

    Returns
    -------
    float
        The redshift corresponding the the DM using the provided method and 
        cosmology.
    """

    table_name = "_".join([method, cosmology]) 
    table_extension = "npy"
    filename = ".".join([table_name, table_extension])
    lookup_table = utils.load_lookup_table(filename)
    z = lookup_table(dm)[()]
    return z



