from __future__ import division, print_function, absolute_import

import os
import numpy as np
from six import PY3, PY2
import scipy.interpolate as interpolate

from fruitbat.methods import avaliable_methods, method_functions


def load(name, data_dir='data'):
    """
    Opens a saved `.npy` file containing an interpolated 1D function.

    Parameters
    ----------
    name : str
        The name of the file to load.

    data_dir : str, optional
        The directory containing the data. The whole path must be
        specified except if :attr:`data_dir` == 'data' then it will
        search in the `data` subdirectory of the source code.
        Default: 'data'

    Returns
    -------
    table: :obj:`scipy.interpolate.interp1d`
        The lookup table containing the interpolated datset.

    Example
    -------
    >>> table = fruitbat.table.load('Zhang2018_Planck18.npy')
    >>> table(1100)
    array(1.22058572)
    """
    if data_dir == 'data':
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
    file = os.path.join(data_dir, name)
    return np.load(file)[()]


def create(filename, method, output_dir='data', zmin=0, zmax=20,
           num_samples=10000, **kwargs):
    """
    Creates an interpolated 1D redshift lookup table which can be read
    in  using :func:`~fruitbat.table.load`.

    Parameters
    ----------
    filename : str
        The name of the output file.

    method : str
        The DM-redshift relation to assume when creating the table.

    output_dir : str
        The path of the output directory. If ``output_dir = 'data'``,
        then created table will created in the same directory with
        the builtin tables and will be found when using functions
        such as `~:meth:fruitbat.Frb.calc_redshift()`.
        Default: 'data'

    zmin : int or float, optional
        The minimum redshift in the table. Default: 0

    zmax : int or float, optional
        The maximum redshift in the table. Default: 20

    num_samples : int, optional
        The number of dispersion measure samples to perform before
        interpolation. Default: 10000

    Keyword Arguments
    -----------------
    cosmo : An instance of :obj:`astropy.cosmology`, optional
        The cosmology to assume when calculating the outputs of the
        table. Required when creating new tables of ``'Ioka2003'``,
        ``'Inoue2004'``, ``'Zhang2018'``.

    free_elec : float or int, optional
        The amount of free electrons per proton mass in the Universe.
        This applies when using ``'Zhang2018'``. Must be between 0
        and 1. Default: 0.875.

    f_igm : float or int, optional
        The fraction of baryons in the intergalactic medium. This
        applies when using ``'Zhang2018'``. Must be between 0 and 1.
        Default: 0.83

    Generates
    ---------
    custom_filename: A binary file in .npy format
        The new lookup table.

    Warning
    -------
    Generating lookup tables is only avaliable when using ``fruitbat``
    in Python 3.
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
                zmin=zmin, zmax=zmax, num_samples=num_samples, **kwargs)

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


def _perform_interpolation(dm_func, zmin=0, zmax=20,
                           num_samples=10000, **kwargs):
    """
    Calculated the dispersion meaasure - redshift relation for a given
    number of samples and returns the interpolated relation.

    Parameters
    ----------
    dm_func : function
        A function that calculates dispersion measure from redshift.
        The first argument of the function must be ``z``. The remaining
        parameters should be keyword only arguments.

    zmin : int or float, optional
        The minimum redshift in the table. Default: 0

    zmax : int or float, optional
        The maximum redshift in the table. Default: 20

    num_samples : int, optional
        The number of dispersion measure samples to perform before
        interpolation. Default: 10000

    Keyword Arguments
    -----------------
    All additional arguments to the dispersion measure function should
    be passed as kewword arguments. The following are example kwargs
    used in teh buildin dispersion measure functions.

    cosmo : An instance of :obj:`astropy.cosmology`, optional
        The cosmology to assume when calculating the outputs of the
        table. Required when creating new tables of ``'Ioka2003'``,
        ``'Inoue2004'``, ``'Zhang2018'``.

    free_elec : float or int, optional
        The amount of free electrons per proton mass in the Universe.
        This applies when using ``'Zhang2018'``. Must be between 0
        and 1. Default: 0.875.

    f_igm : float or int, optional
        The fraction of baryons in the intergalactic medium. This
        applies when using ``'Zhang2018'``. Must be between 0 and 1.
        Default: 0.83

    Returns
    -------
    interp: An instance method of :obj:`scipy.interpolate.interp1d`
        The interpolated lookup table.
    """
    z_vals = np.linspace(zmin, zmax, num_samples)
    dm_vals = np.array([dm_func(z, **kwargs) for z in z_vals])
    interp = interpolate.interp1d(dm_vals, z_vals)

    return interp
