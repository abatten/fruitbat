import numpy as np

__all__ = ['dm_to_redshift']


def dm_to_redshift(dm, dm_err=None):

    z = np.sqrt(dm) - 1

    if dm_err:
        z_err = np.sqrt(dm_err) - 1
    else:
        z_err = 0.0

    return z, z_err


def calc_redshift(dm, dm_err=None, method='batten2019'):
    """
    Returns the redshift of a given dispersion measure using a specified method.

    Parameters
    ----------
    dm : float
        Dispersion Measure.

    dm_err : float or None
        The uncertianty in the dispersion measure. If dm_err is `None`, 

    method : string
        The z-DM relation to use to calculate the redshift. Avaliable methods 
        are: `batten2019`

    Returns
    -------
    z : float
        Redshift

    z_err : float
        The uncertianty in the redshift estimation. If dm_err is `None` then
        `z_err` = 0.
    """
   
    valid_methods = ['batten2019']

    if method not in valid_methods:
        raise ValueError("""Method '{0}' is not a valid method for calc_redshift. 
                            Valid methods are: 'batten2019'""".format(method))

    if method == 'batten2019':
        z, z_err = _batten2019(dm, dm_err)

    return z, z_err

def _batten2019(dm, dm_err=None):
    return 12.0, 4.0
