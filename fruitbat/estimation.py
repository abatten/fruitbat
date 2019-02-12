"""
Estimation
==========

The estimation module
"""

import numpy as np
import os

# __all__ = ['dm_to_redshift']


def calc_redshift(dm, dm_uncert=None, method='batten2019', 
                  cosmology='planck2018+bao'):
    """
    Returns the redshift of a given dispersion measure using a
    specified DM-z relation.

    Parameters
    ----------
    dm : float
        Dispersion Measure.

    dm_uncert : float or None
        The uncertianty in the dispersion measure. If ``dm_err`` is `None`,

    method : string, optional
        The z-DM relation to use to calculate the redshift.
        Avaliable methods are:
        ``'batten2019'``, ``'inoue2004'``, ``'ioka2003'``
    
    cosmology : string, optional
        Avaliable cosmologies:

    Returns
    -------
    z : float
        Redshift

    z_err : float
        The uncertianty in the redshift estimation. If `dm_uncert` is `None`
        then `z_err` = 0.
    """

    valid_methods = ['batten2019', 'zhang2018', 'inoue2004', 'ioka2003']

    if method not in valid_methods:
        raise ValueError("""Method '{m}' is not a valid method.
            Valid methods are: {valid}""".format(m=method, valid=valid_methods))

    if method == 'batten2019':
        z, z_uncert = _dm_to_z_batten2019(dm, dm_uncert)

    elif method == 'inoue2004':
        z, z_uncert = _dm_to_z_inoue2004(dm, dm_uncert, cosmology)

    elif method == 'ioka2003':
        z, z_uncert = _dm_to_z_ioka2003(dm, dm_uncert)

    return z, z_uncert


def _dm_to_z_batten2019(dm, dm_err=None):
    """
    Calculates a redshft from a dispersion measure using the DM-z
    relation from Batten, A. J. 2019, ....
    """
    return 12.0, 4.0


def _dm_to_z_inoue2004(dm, dm_uncert=None, cosmology="planck2018+bao"):
    """
    Calculates a redshift from a dispersion measure using the DM-z
    relation from Inoue, S. 2004, MNRAS, 348(3), 999–1008.

    Pararameters
    ------------

    Returns
    -------
    z : float
        The redshift of the FRB.
    z_uncert : tuple of floats
        The uncertainty in the redshift of the FRB.
    """

    cosmo_dict = {
        "wmap2012": "inoue2004_wmap2012.npy",
        "planck2015": "inoue2004_planck2015.npy",
        "planck2018": "inoue2004_planck2018.npy", 
        "planck2018+bao": "inoue2004_planck2018_bao.npy",
    }

    lookup_table = load_lookup_table(cosmo_dict[cosmology])

    z, dz_low, dz_high = _get_redshift_from_table(lookup_table, dm, dm_uncert) 

    return z, (dz_low, dz_high)


def _dm_to_z_ioka2003(dm, dm_err=None):
    """
    Calculates a redshift from a dispersion measure using the DM-z
    relation from Ioka, K. 2003, ApJL, 598, L79
    """
    return 100.0, 5.0


def load_lookup_table(filename, data_dir=''):
    """
    Parameters
    ----------
    """
    if data_dir == '':
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
    interp_file = os.path.join(data_dir, filename) 
    return np.load(interp_file)[()]


def _get_redshift_from_table(table, dm, d_dm):

    z = table(dm)[()]

    dz_low = abs(z - table(dm - d_dm)[()])
    dz_high = abs(z - table(dm + d_dm)[()])

    return z, dz_low, dz_high
