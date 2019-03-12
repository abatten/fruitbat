"""
Utility functions for Fruitbat
"""

import os

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import astropy.constants as const
import astropy.units as u

__all__ = ["create_lookup_table", "load_lookup_table"]


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
    :obj:`scipy.interpolate.interp1d`
        Function
    """
    if data_dir == 'data':
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
    file = os.path.join(data_dir, filename)
    return np.load(file)[()]


def create_lookup_table(filename, method, cosmo, zmin=0, zmax=20,
                        num_samples=1e4, **kwargs):
    """
    Creates a lookup table

    Parameters
    ----------
    filename : str
        The name

    method : str

    cosmo : dict
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
        This only applies when using ``'zhang2018'``. Must be between 0 and 1.
        Default: 0.875.

    f_igm : float or int, optional
        The fraction of baryons in the intergalactic medium. This only
        applies when using ``'zhang2018'``. Must be between 0 and 1. 
        Default: 0.83

    Generates
    ---------
    filename.npy:

    """

#    method_table_dict = {"inoue2004": _create_lookup_table_inoue2004,
#                         "zhang2018": _create_lookup_table_zhang2018,
#                         "ioka2003": _create_lookup_table_ioka2003
#                         }

    dm_method_dict = {
        "ioka2003": _calc_dm_ioka2003,
        "inoue2004": _calc_dm_inoue2004,
        "zhang2018": _calc_dm_zhang2018
    }

    interp = _perform_interpolation(dm_func=dm_method_dict[method],
                                    cosmo=cosmo,
                                    zmin=zmin, zmax=zmax,
                                    num_samples=num_samples,
                                    **kwargs)

    _save_lookup_table(interp, filename)


#    method_table_dict[method](filename, cosmology, zmin=zmin, zmax=zmax,
#                              num_samples=num_samples)


def _fz_integrand(z, cosmo):
    """
    Calculate the integrand for a given redshift and cosmology. This integral
    appears in zhang2018, inoue2004 and ioka2003.

    Parameters
    ----------
    z : float
        The redshift to evaluate the integral

    cosmology
    """

    w = cosmo.w(z)

    top = 1 + z
    bot = cosmo.Om0 * (1 + z)**3 + cosmo.Ode0 * (1 + z)**(3 + 3 * w)
    return top / np.sqrt(bot)


def _calc_dm_ioka2003(z, cosmo, zmin=0):
    """
    Calculate the dispersion measure from a redshift given a cosmology
    using the Ioka (2003) relation.

    Parameters
    ----------
    z: float or int
        The redshift
    cosmo:

    Returns
    -------
    dm : float
        The dispersion measure at the redhsift.
    """
    # Calculate Ioka 2003 DM coefficient
    coeff_top = 3 * const.c * cosmo.H0 * cosmo.Ob0
    coeff_bot = 8 * np.pi * const.G * const.m_p
    coeff = coeff_top / coeff_bot

    coeff = coeff.to("pc cm**-3")

    dm = coeff * integrate.quad(_fz_integrand, zmin, z, args=(cosmo))[0]

    return dm.value


def _calc_dm_inoue2004(z, cosmo, zmin=0):
    """
    Calculate the dispersion measure from a redshift given a cosmology
    using the Inoue (2004) relation.

    Parameters
    ----------
    z :
    cosmology :
    zmin: float or int, optional
        The minimum redshift

    Returns
    -------
    dm : float
        The dispersion measure at the redshift.

    """
    inoue_n_e_0 = 9.2e-10 * ((u.Mpc**2 * u.s**2) / (u.km**2 * u.cm**3))

    coeff = inoue_n_e_0 * const.c * cosmo.Ob0 * cosmo.H0
    coeff = coeff.to("pc cm**-3")

    dm = coeff * integrate.quad(_fz_integrand, zmin, z, args=(cosmo))[0]

    return dm.value


def _calc_dm_zhang2018(z, cosmo, zmin=0, **kwargs):
    """
    Calculate the dispersion measure from a redshift given a cosmology
    using the Zhang (2018) relation.

    Parameters
    ----------
    z :
    cosmology :
    zmin: float or int, optional
        The minimum redshift

    Returns
    -------
    dm : float
        The dispersion measure at the redshift.
    """

    if "f_igm" in kwargs:
        if kwargs["f_igm"] <= 1 and kwargs["f_igm"] >= 0:
            f_igm = kwargs["f_igm"]
        else:
            raise ValueError("f_igm must be between 0 and 1.")
    else:
        f_igm = 0.83

    if "free_elec" in kwargs:
        if kwargs["free_elec"] <= 1 and kwargs["free_elec"] >= 0:
            free_elec = kwargs["free_elec"]
        else:
            raise ValueError("free_elec must be between 0 and 1.")
    else:
        free_elec = 0.875

    coeff_top = 3 * const.c * cosmo.H0 * cosmo.Ob0
    coeff_bot = 8 * np.pi * const.G * const.m_p
    coeff = coeff_top / coeff_bot

    coeff = coeff.to("pc cm**-3")
    dm = coeff * f_igm * free_elec * integrate.quad(_fz_integrand, zmin, z,
                                                    args=(cosmo))[0]
    return dm.value

#def _create_lookup_table_ioka2003(filename, cosmo, zmin=0, zmax=20,
#                                  num_samples=1e4):
#    """
#    Creates an interpolated 1D DM-z look up table using the Ioka (2003)
#    relation and a given cosmology.
#
#    Parameters
#    ----------
#    filename: str
#
#    cosmo: dict
#
#    zmin: int or float, optional
#        The minimum redshift for the table. The output table will
#        not be able to estimate redshifts lower than this value. Default: 0
#
#    zmax: int or float, optional
#        The maximum redshift for the table. The output table will
#        not be able to estimate redshifts higher than this value. Default: 20
#
#    num_samples: int, optional
#        Default: 10000
#
#    Returns
#    -------
#    None
#    """
#
#    def _calc_dm(z, cosmo):
#        """
#        Calculate the dispersion measure from a redshift given a cosmology
#        using the Ioka (2003) relation.
#        """
#        # Check that the user has provided all the required values
#        #cosmo_required_keys = ["HO", "Omega_m", "Omega_b", "Omega_L"]
#        #_check_keys_in_dict(cosmo, cosmo_required_keys)
#
#        coeff_top = 3 * const.c * cosmo.H0 * cosmo.Ob0
#        coeff_bot = 8 * np.pi * const.G * const.m_p
#        coeff = coeff_top / coeff_bot
#        coeff = coeff.to("pc cm**-3")
#
#        dm = coeff * integrate.quad(_fz_integrand, 0, z, args=(cosmo))[0]
#
#        return dm.value
#
#    interp = _perform_interpolation(dm_func=_calc_dm, cosmology=cosmo,
#                                    zmin=zmin, zmax=zmax,
#                                    num_samples=num_samples)
#
#    _save_lookup_table(interp, filename)
#
#
#
#
#
#def _create_lookup_table_inoue2004(filename, cosmo, zmin=0, zmax=20,
#                                   num_samples=1e4):
#    """
#    Creates an interpolated 1D DM-z look up table using the Inoue (2004)
#    relation and a given cosmology.
#
#    Parameters
#    ----------
#    filename: str
#
#    cosmo: dict
#
#    zmin: int or float, optional
#        The minimum redshift for the table. The output table will
#        not be able to estimate redshifts lower than this value. Default: 0
#
#    zmax: int or float, optional
#        The maximum redshift for the table. The output table will
#        not be able to estimate redshifts higher than this value. Default: 20
#
#    num_samples: int, optional
#        Default: 10000
#
#    Returns
#    -------
#    None
#    """
#
#    def _calc_dm(z, cosmo):
#        """
#        Calculate the dispersion measure from a redshift given a cosmology
#        using the Inoue (2004) relation.
#        """
#        # Check that the user has provided all the required values
#        #cosmo_required_keys = ["HO", "Omega_m", "Omega_b", "Omega_L"]
#        #_check_keys_in_dict(cosmo, cosmo_required_keys)
#
#        inoue_n_e_0 = 9.2e-10 * ((u.Mpc**2 * u.s**2) / (u.km**2 * u.cm**3))
#
#        coeff = inoue_n_e_0 * const.c * cosmo.Ob0 * cosmo.H0
#        coeff = coeff.to("pc cm**-3")
#
#        dm = coeff * integrate.quad(_fz_integrand, 0, z, args=(cosmo))[0]
#
#        return dm.value
#
#    interp = _perform_interpolation(dm_func=_calc_dm, cosmology=cosmo,
#                                    zmin=zmin, zmax=zmax,
#                                    num_samples=num_samples)
#
#    _save_lookup_table(interp, filename)
#
#
#def _create_lookup_table_zhang2018(filename, cosmo, zmin=0, zmax=20,
#                                   num_samples=1e4, f_igm=0.83,
#                                   free_elec=0.875):
#    """
#    Creates an interpolated 1D DM-z look up table using the Zhang (2018)
#    relation and a given cosmology.
#
#    Parameters
#    ----------
#    filename: str
#
#    cosmo: dict
#
#    zmin: int or float, optional
#        Default: 0
#
#    zmax: int or float, optional
#        Default: 30
#
#    num_samples: int, optional
#        Default: 100000
#
#    f_igm: float, optional
#        Default: 0.83
#
#    free_elec: float, optional
#        Default: 0.875
#    """
#
#    def _calc_dm(z, cosmo, approx=False, *args, **kwargs):
#
#        # Check that the user has provided all the required values
#        #cosmo_required_keys = ["HO", "Omega_m", "Omega_b", "Omega_L"]y
#        #_check_keys_in_dict(cosmo, cosmo_required_keys)
#
#        coeff_top = 3 * const.c * cosmo.H0 * cosmo.Ob0
#        coeff_bot = 8 * np.pi * const.G * const.m_p
#        coeff = coeff_top / coeff_bot
#
#        coeff = coeff.to("pc cm**-3")
#        dm = coeff * f_igm * free_elec * integrate.quad(_fz_integrand, 0, z,
#                                                        args=(cosmo))[0]
#        return dm.value
#
#    interp = _perform_interpolation(dm_func=_calc_dm, cosmology=cosmo,
#                                    zmin=zmin, zmax=zmax,
#                                    num_samples=num_samples)
#
#    _save_lookup_table(interp, filename)


def _perform_interpolation(dm_func=None, zmin=None, zmax=None, cosmo=None,
                           num_samples=None, **kwargs):
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


def _check_keys_in_dict(dictionary, keys):
    """
    Checks that a list of keys exist in a dictionary.

    Parameters
    ----------
    dictionary: dict


    keys: list of strings
        The keys that the dictionary must have.

    Returns
    -------
    bool:
        True
    """
    if not all(key in dictionary for key in keys):
        raise KeyError("Dictionary missing key values."
                       "Requires: {}".format(keys))
    return True


def _save_lookup_table(table, filename):
    """
    Saves the lookup table to a .npy file
    """
    np.save(filename, table)
    return None
