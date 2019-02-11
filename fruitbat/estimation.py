import numpy as np

# __all__ = ['dm_to_redshift']


def calc_redshift(dm, dm_err=None, method='batten2019'):
    """
    Returns the redshift of a given dispersion measure using a
    specified DM-z relation.

    Parameters
    ----------
    dm : float
        Dispersion Measure.

    dm_err : float or None
        The uncertianty in the dispersion measure. If ``dm_err`` is `None`,

    method : string
        The z-DM relation to use to calculate the redshift.
        Avaliable methods are:
        ``'batten2019'``, ``'inoue2004'``, ``'ioka2003'``

    Returns
    -------
    z : float
        Redshift

    z_err : float
        The uncertianty in the redshift estimation. If dm_err is `None`
        then ``z_err`` = 0.
    """

    valid_methods = ['batten2019', 'inoue2004', 'ioka2003']

    if method not in valid_methods:
        raise ValueError("""Method '{m}' is not a valid method.
            Valid methods are: {valid}""".format(m=method, valid=valid_methods))

    if method == 'batten2019':
        z, z_err = _dm_to_z_batten2019(dm, dm_err)

    elif method == 'inoue2004':
        z, z_err = _dm_to_z_inoue2004(dm, dm_err)

    elif method == 'ioka2003':
        z, z_err = _dm_to_z_ioka2003(dm, dm_err)

    return z, z_err


def _dm_to_z_batten2019(dm, dm_err=None):
    """
    Calculates a redshft from a dispersion measure using the DM-z
    relation from Batten, A. J. 2019, ....
    """
    return 12.0, 4.0


def _dm_to_z_inoue2004(dm, dm_err=None,
                       cosmo={'HO': 70, 'Omega_m': 0.3, 'Omega_L': 0.7}):
    """
    Calculates a redshift from a dispersion measure using the DM-z
    relation from Inoue, S. 2004, MNRAS, 348(3), 999–1008.
    """
    return 50.0, 2.0


def _dm_to_z_ioka2003(dm, dm_err=None):
    """
    Calculates a redshift from a dispersion measure using the DM-z
    relation from Ioka, K. 2003, ApJL, 598, L79
    """
    return 100.0, 5.0
