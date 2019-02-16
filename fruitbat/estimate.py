from . import utils

# __all__ = ['dm_to_redshift']
from ._fruitbatstrings import (docstr_sub, _cosmo_doc,
                               _methods_doc, _dm_units_doc)

@docstr_sub(dmunits=_dm_units_doc, methods=_methods_doc, cosmo=_cosmo_doc)
def redshift(dm, dm_uncert=0.0, method='inoue2004', cosmology='planck2018'):
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

    valid_methods = ['batten2019', 'zhang2018', 'inoue2004', 'ioka2003']

    if method not in valid_methods:
        raise ValueError("""Method '{m}' is not a valid method.
            Valid methods are: %(cosmos)s""")

    if method == 'batten2019':
        z = _redshift_batten2019(dm, dm_uncert)

    elif method == 'inoue2004':
        z = _redshift_inoue2004(dm, dm_uncert, cosmology)

    elif method == 'ioka2003':
        z = _redshift_ioka2003(dm, dm_uncert)

    elif method == 'zhang2018':
        z = _redshift_zhang2018(dm, dm_uncert, cosmology)

    return z


def _redshift_batten2019(dm, dm_uncert=0.0):
    """
    Calculates a redshft from a dispersion measure using the DM-z
    relation from Batten, A. J. 2019, ....
    """
    return 12.0

@docstr_sub(dm_units=_dm_units_doc)
def _redshift_inoue2004(dm, dm_uncert=0.0, cosmology="planck2018"):
    """
    Calculates a redshift from a dispersion measure using the DM-z
    relation from Inoue2004_

    Parameters
    ----------
    dm : float
        Units: %(dm_units)s

    dm_uncert : float, optional
        Default: :0.0

    cosmology : str, optional
        Default: `'planck2018'`

    Returns
    -------
    z : float
        The redshift of the FRB.
    z_uncert : tuple of floats
        The uncertainty in the redshift of the FRB.

    Notes
    -----
    Inoue2004_ presents the following dispersion measure and redshift relation.
   
    .. math ::
        \\mathrm{DM}(z) =0.92\\times 10^{-5} \\Omega_b h^2 c \\int_0^z
        \\frac{1 + z'}{H_0\\left(\\Omega_m (1 + z)^3 +
        \\Omega_\\Lambda\\right)^{1/2}} dz'


    .. _Inoue2004: http://adsabs.harvard.edu/abs/2004MNRAS.348..999I
    """

    cosmo_dict = {
        "wmap2013": "inoue2004_wmap2013.npy",
        "planck2015": "inoue2004_planck2015.npy",
        "planck2018": "inoue2004_planck2018.npy", 
        "planck2018+bao": "inoue2004_planck2018_bao.npy",
    }

    lookup_table = utils.load_lookup_table(cosmo_dict[cosmology])

    z = _get_redshift_from_table(lookup_table, dm, dm_uncert) 
    #z, dz_low, dz_high = _get_redshift_from_table(lookup_table, dm, dm_uncert) 
    return z#, (dz_low, dz_high)


def _redshift_ioka2003(dm, dm_uncert=0.0):
    """
    Calculates a redshift from a dispersion measure using the DM-z
    relation from Ioka, K. 2003, ApJL, 598, L79
    """
    return 100.0


def _redshift_zhang2018(dm, dm_uncert=0.0, cosmology="planck2018"):
    """
    Calculates a redshift from a dispersion measure using the DM-z
    relation from Zhang2018_


    .. _Zhang2018: https://ui.adsabs.harvard.edu/#abs/arXiv:1808.05277
    """
    return 187.0


def _get_redshift_from_table(table, dm, dm_uncert):

    z = table(dm)[()]
    #dz_low = abs(z - table(dm - dm_uncert)[()])
    #dz_high = abs(z - table(dm + dm_uncert)[()])

    return z#, dz_low, dz_high
