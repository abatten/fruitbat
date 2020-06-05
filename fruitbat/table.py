"""
table
=====

"""
import os
import numpy as np
import scipy.interpolate as interpolate
import h5py

from fruitbat import utils
from fruitbat.methods import available_methods, method_functions


__all__ = ["create", "get_z_from_table", "get_table_path"]


def create(method, output_dir='data', filename=None, zmin=0, zmax=20,
           num_samples=10000, **method_kwargs):
    """
    Creates lookup table


     which can be read
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


    Returns
    -------
    string
        The path to the generated hdf5 file containing the table data.

    Generates
    ---------
    A hdf5 file containing datasets for `'DM'` and 'z'`.


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

    # The filename that the user provides can not be one of these as it would
    # over write the built in datasets.
    restricted_filenames = set([
        utils.get_path_to_file_from_here("Ioka2003.hdf5", subdirs=["data"]),
        utils.get_path_to_file_from_here("Inoue2004.hdf5", subdirs=["data"]),
        utils.get_path_to_file_from_here("Zhang2018.hdf5", subdirs=["data"]),
        utils.get_path_to_file_from_here("Batten2020.hdf5", subdirs=["data"]),
    ])

    # Generate the output filename.
    if (filename is not None) and (filename not in restricted_filenames):
        filename = "{}.hdf5".format(filename)
    else:
        filename = "{}.hdf5".format(method)

    # Generate the DM and z values for the method
    dm_function = method_functions()[method]
    z_vals = np.linspace(zmin, zmax, num_samples)
    dm_vals = np.array([dm_function(z, **method_kwargs) for z in z_vals])

    if output_dir == "data":
        output_file = utils.get_path_to_file_from_here(filename, subdirs=["data"])

    else:
        output_file = os.path.join(output_dir, filename)


    with h5py.File(output_file, "w") as new_table:
        new_table.create_group("Header")

        if "cosmo" in method_kwargs:
            cosmo = method_kwargs["cosmo"]
            add_cosmo_params_to_dataset_attrs = True

            # Use a default cosmology name if the prodived one doesnt have a name
            if not cosmo.name:
                cosmo.name = "cosmology"

            dataset_group_name = cosmo.name
        else:
            add_cosmo_params_to_dataset_attrs = False
            dataset_group_name = method

        new_table.create_group(dataset_group_name)

        new_table[dataset_group_name].create_dataset("DM", data=dm_vals, dtype=np.float)
        new_table[dataset_group_name]["DM"].attrs["Units"] = "pc cm**-3"
        new_table[dataset_group_name]["DM"].attrs["VarDesc"] = "Dispersion Measure"

        new_table[dataset_group_name].create_dataset("z", data=z_vals, dtype=np.float)
        new_table[dataset_group_name]["z"].attrs["Units"] = "Dimensionless"
        new_table[dataset_group_name]["z"].attrs["VarDesc"] = "Redshift"

        if add_cosmo_params_to_dataset_attrs:
            new_table[dataset_group_name].attrs["Flat"] = True
            new_table[dataset_group_name].attrs["H0"] = cosmo.H0.value
            new_table[dataset_group_name].attrs["OmegaBaryon0"] = cosmo.Ob0
            new_table[dataset_group_name].attrs["OmegaLambda0"] = cosmo.Ode0
            new_table[dataset_group_name].attrs["OmegaMatter0"] = cosmo.Om0
            new_table[dataset_group_name].attrs["t0"] = cosmo.age(0).value

    return output_file









# def load(name, data_dir='data'):
#     """
#     Opens a saved `.npz` file containing 'dm' and 'z' arrays.

#     Parameters
#     ----------
#     name : str
#         The name of the file to load.

#     data_dir : str, optional
#         The directory containing the data. The whole path must be
#         specified except if :attr:`data_dir` == 'data' then it will
#         search in the `data` subdirectory of the source code.
#         Default: 'data'

#     Returns
#     -------
#     table: :obj:`numpy.lib.npyio.NpzFile`
#         The lookup table containing the 'dm' and 'z' arrays.

#     Example
#     -------
#     >>> table = fruitbat.table.load('Zhang2018_Planck18.npz')
#     >>> table["dm"]
#     array([0.00000000e+00, 1.62251609e+00, 3.24675204e+00, ...,
#            1.00004587e+04, 1.00010926e+04, 1.00017266e+04])
#     >>> table["z"]
#     array([0.00000000e+00, 2.00020002e-03, 4.00040004e-03, ...,
#            1.99959996e+01, 1.99979998e+01, 2.00000000e+01])

#     """
#     if data_dir == 'data':
#         data_dir = os.path.join(os.path.dirname(__file__), 'data')

#     filename = os.path.join(data_dir, name)
#     return np.load(filename)


def get_table_path(filename, datadir="data"):
    """

    Parameters
    ----------
    filename:

    datadir:

    Returns
    -------
    path: string
        The path to the data file
    """

    if filename in available_methods():
        filename += ".hdf5"

    if datadir == "data":
        path = utils.get_path_to_file_from_here(filename, subdirs=["data"])

    return path


def get_z_from_table(dm, table, cosmology=None):
    """
    Calculates the redshift from a dispersion measure by interpolating
    a lookup table

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
    with h5py.File(table, "r") as data:
        if cosmology:
            dataset_group_name = cosmology
        else:
            dataset_group_name = os.path.basename(table).split(".")[0]
        interp = interpolate.interp1d(data[dataset_group_name]["DM"], data[dataset_group_name]["z"])
    return interp(dm)
