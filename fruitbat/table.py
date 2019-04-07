from __future__ import division, print_function, absolute_import

import os
import numpy as np
from six import PY3, PY2
import scipy.interpolate as interpolate

from fruitbat.methods import avaliable_methods, method_functions


def load(filename, data_dir='data'):
    """
    Opens a saved `.npy` file containing an interpolated 1D function.

    Parameters
    ----------
    filename : str
        The name of the file to load.

    data_dir : str, optional
        The directory containing the data. The whole path must be specified
        except if :attr:`data_dir` == 'data' then it will search in the
        `data` subdirectory of the source code. Default: 'data'

    Returns
    -------
    :obj:`scipy.interpolate.interp1d`
        Function
    """
    if data_dir == 'data':
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
    file = os.path.join(data_dir, filename)
    return np.load(file)[()]


def create(filename, method, cosmo, output_dir='data', zmin=0, zmax=20, 
           num_samples=10000, **kwargs):
    """
    Creates an interpolated 1D redshift lookup table which can be read in
    using :func:`load_lookup_table`.

    Parameters
    ----------
    filename : str
        The name of the output file.

    method : str
        The DM-z relation to assume when creating the table.

    cosmo : :obj:`astropy.cosmology` or None
        A dictionary containing the cosmology parameters for the Hubble
        constant, Matter density, Baryon density and Dark Energy density.
        `cosmology` must contain values for the following keys:
        ``'H0'``, ``'Omega_m'``, ``'Omega_b'``, ``'Omega_L'``

    zmin : int or float, optional
        The minimum redshift in the table. Default: 0

    zmax : int or float, optional
        The maximum redshift in the table. Default: 20

    num_samples : int, optional
        The number of dispersion measure samples to perform before
        interpolation. Default: 10000

    Keyword Arguments
    -----------------
    free_elec : float or int, optional
        The amount of free electrons per proton mass in the Universe.
        This only applies when using ``'Zhang2018'``. Must be between 0 and 1.
        Default: 0.875.

    f_igm : float or int, optional
        The fraction of baryons in the intergalactic medium. This only
        applies when using ``'Zhang2018'``. Must be between 0 and 1.
        Default: 0.83

    Generates
    ---------
    filename.npy:

    Warning
    -------
    Generating lookup tables is only avaliable when using **fruitbat** in
    Python 3.

    """
    if PY3:

        if method not in avaliable_methods():
            err_msg = ("{} is not a valid method."
                       "The currently defined methods "
                       "are: {}".format(method, avaliable_methods()))
            raise ValueError(err_msg)

        else:
            function = method_functions()[method]
            lookup_table = _perform_interpolation(function,
                                                  cosmo=cosmo,
                                                  zmin=zmin, zmax=zmax,
                                                  num_samples=num_samples,
                                                  **kwargs)

            output_name = "custom_{}.npy".format(filename)
            if output_dir == 'data':
                output_file = os.path.join(os.path.dirname(__file__),
                                           'data', output_name)
            else:
                output_file = os.path.join(output_dir, output_name)

            np.save(output_file, lookup_table)

    elif PY2:
        raise SystemError("""table.create is a python 3 only feature.
                          A lookup table is in reality a pickled
                          instance method of scipy.interpolate.interp1D.
                          Since pickling an instance method in Python 2
                          is not supported, this is restricted to
                          Python 3 only.""")


def _perform_interpolation(dm_func, cosmo=None, zmin=0, zmax=20,
                           num_samples=10000, **kwargs):
    """
    Parameters
    ----------

    Returns
    -------


    """
    z_vals = np.linspace(zmin, zmax, num_samples)
    dm_vals = np.array([dm_func(z, cosmo, **kwargs) for z in z_vals])
    interp = interpolate.interp1d(dm_vals, z_vals)

    return interp
