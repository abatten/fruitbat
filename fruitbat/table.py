from __future__ import division, print_function, absolute_import

import os
import numpy as np
import scipy.interpolate as interpolate

from fruitbat.methods import available_methods, method_functions

__all__ = ["create", "load", "get_z_from_table"]


def create(method, output_dir='data', filename=None, zmin=0, zmax=20,
           num_samples=10000, **kwargs):
    """
    Creates an interpolated 1D redshift lookup table which can be read
    in  using :func:`~fruitbat.table.load`.

    Parameters
    ----------
    method : str
        The DM-redshift relation to assume when creating the table.

    output_dir : str, optional
        The path of the output directory. If ``output_dir = 'data'``,
        then created table will created in the same directory with
        the builtin tables and will be found when using functions
        such as :meth:`~fruitbat._frb.calc_redshift()`.
        Default: 'data'

    filename : str, optional
        The output filename. If ``name=None`` then the filename will
        become custom_method. Default: *None*

    zmin : int or float, optional
        The minimum redshift in the table. Default: 0

    zmax : int or float, optional
        The maximum redshift in the table. Default: 20

    num_samples : int, optional
        The number of sample dispersion measure and redshifts.
        Default: 10000

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
    A compressed file in .npz format containing the arrays ``'dm'``
    and ``'z'``.

    Example
    -------
    >>> def simple_dm(z):
        dm = 1200 * z
        return dm
    >>> fruitbat.add_method("simple_dm", simple_dm)
    >>> fruitbat.table.create("simple_dm")
    >>> frb = fruitbat.Frb(1200)
    >>> frb.calc_redshift(method="simple_dm")
    <Quantity 1.>

    """

    if method not in available_methods():
        err_msg = ("{} is not a valid method."
                   "The currently defined methods "
                   "are: {}".format(method, available_methods()))
        raise ValueError(err_msg)

    else:
        dm_function = method_functions()[method]

        z_vals = np.linspace(zmin, zmax, num_samples)
        dm_vals = np.array([dm_function(z, **kwargs) for z in z_vals])

        if filename is not None:
            output_name = filename
        else:
            output_name = "custom_{}".format(method)

        if output_dir == 'data':
            output_file = os.path.join(os.path.dirname(__file__),
                                       'data', output_name)
        else:
            output_file = os.path.join(output_dir, output_name)

        np.savez(output_file, dm=dm_vals, z=z_vals)


def load(name, data_dir='data'):
    """
    Opens a saved `.npz` file containing 'dm' and 'z' arrays.

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
    table: :obj:`numpy.lib.npyio.NpzFile`
        The lookup table containing the 'dm' and 'z' arrays.

    Example
    -------
    >>> table = fruitbat.table.load('Zhang2018_Planck18.npz')
    >>> table["dm"]
    array([0.00000000e+00, 1.62251609e+00, 3.24675204e+00, ...,
           1.00004587e+04, 1.00010926e+04, 1.00017266e+04])
    >>> table["z"]
    array([0.00000000e+00, 2.00020002e-03, 4.00040004e-03, ...,
           1.99959996e+01, 1.99979998e+01, 2.00000000e+01])

    """
    if data_dir == 'data':
        data_dir = os.path.join(os.path.dirname(__file__), 'data')

    filename = os.path.join(data_dir, name)
    return np.load(filename)


def get_z_from_table(dm, table):
    """
    Calculates the redshift from a dispersion meausre by interpolating
    a lookup table.

    Parameters
    ----------
    dm: float
        The input dispersion measure

    table: :obj:`numpy.lib.npyio.NpzFile`
        The lookup table with ``'dm'`` and ``'z'`` arrays.

    Returns
    -------
    z: float
        The redshift corresponding to the input disperison measure.

    Example
    -------
    >>> table = fruitbat.table.load('Zhang2018_Planck18.npz')
    >>> fruitbat.get_z_from_table(1000, table)
    1.1087964578507539

    """
    interp = interpolate.interp1d(table["dm"], table["z"])
    return interp(dm)
