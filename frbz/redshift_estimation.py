import numpy as np

__all__ = ['dm_to_redshift']


def dm_to_redshift(dm, dm_err=None):

    z = np.sqrt(dm) - 1

    if dm_err:
        z_err = np.sqrt(dm_err) - 1
    else:
        z_err = 0.0

    return z, z_err

