"""
Utility functions for Fruitbat
"""

import os

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate


def load_lookup_table(filename, data_dir='data'):
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
    scipy.interpolate.interpolate.interp1d
        Function
    """
    if data_dir == 'data':
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
    file = os.path.join(data_dir, filename)
    return np.load(file)[()]


def create_lookup_table(filename, method, cosmology, zmin=0, zmax=30,
                        num_samples=1e5):
    """
    Creates a lookup table

    Parameters
    ----------
    filename : str

    method : str

    cosmology : dict
        A dictionary containing the cosmology parameters for the Hubble 
        constant, Matter density, Baryon density and Dark Energy density. 
        `cosmology` must contain values for the following keys: 
        ``'H0'``, ``'Omegam'``, ``'Omegab'``, ``'OmegaL'``

    zmin : int or float, optional
        The minimum redshift in the table. Default: 0

    zmax : int or float, optional
        The maximum redshift in the table. Default: 30

    num_samples : int, optional
        The number of dispersion measure samples to perform before 
        interpolation. Default: 100000

    Returns
    -------
    path : str
        The path of the created lookup table

    """

    
    

def _create_lookup_table_inoue2004(filename, cosmo, zmin=0, zmax=30, 
                                   num_samples=1e5):
    """
    Creates an interpolated 1D DM-z look up table using the Inoue (2004)
    relation and a given cosmology.
    """

    def _integrand(z, cosmo):
        """
        Calculate the integrand of the Inoue2004 function
        """
        top = 1 + z
        bot = cosmo["HO"] * (cosmo["Omegam"] * (1 + z)**3 + cosmo["OmegaL"])

        return top / np.sqrt(bot)
        
    def _calc_dm(z, cosmo):
        """
        Calculate the dispersion measure from a redshift given a cosmology
        using the Inoue (2004) relation.
        """
        coeff = 0.92e-5 * cosmo["Omegab"] * (cosmo["HO"]/100)**2
        c = 29979245800  # cm / s

        dm = coeff * c * integrate.quad(_integrand, 0, z, args=(cosmo))[0]

        return dm

    z_vals = np.linspace(zmin, zmax, num_samples)
    dm_vals = np.array([_calc_dm(zi, cosmo) for zi in z_vals])

    interp = interpolate.interp1d(dm_vals, z_vals)

    _save_lookup_table(interp, filename)

    return None


def _save_lookup_table(table, filename):
    np.save(filename, table)
    return None
