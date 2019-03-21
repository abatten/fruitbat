from e13tools import docstring_substitute

from fruitbat import utils
from fruitbat._fruitbatstrings import dm_units_doc
from fruitbat.cosmology import keys as cosmo_keys, builtin

__all__ = ["redshift", "methods"]


def methods(string=False):
    """
    Defines the list of avaliable method keywords.

    Methods currently avaliable: Ioka2003, Inoue2004, Zhang2018

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
    methods = ["Ioka2003", "Inoue2004", "Zhang2018"]

    if string:
        methods = ", ".join(methods)

    return methods


@docstring_substitute(dmunits=dm_units_doc, methods=methods(string=True), 
                      cosmo=cosmo_keys())
def redshift(dm, dm_uncert=0.0, method='Inoue2004', cosmology='Planck18'):
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
        Avaliable methods are: %(methods)s. Default: `'Inoue2004'`

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



