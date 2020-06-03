"""
The collection of utility functions for Fruitbat.
"""
from __future__ import print_function, absolute_import, division
import os
import numpy as np
from scipy import interpolate

__all__ = ["check_keys_in_dict"]


def check_keys_in_dict(dictionary, keys):
    """
    Checks that a list of keys exist in a dictionary.

    Parameters
    ----------
    dictionary: dict
        The input dictionary.

    keys: list of strings
        The keys that the dictionary must contain.

    Returns
    -------
    bool:
        Returns *True* is all required keys exist in the dictionary.
        Otherwise a KeyError is raised.
    """
    if not all(key in dictionary for key in keys):
        raise KeyError("Dictionary missing key values."
                       "Requires: {}".format(keys))
    return True


def get_path_to_file_from_here(filename, subdirs=None):
    """
    Returns the whole path to a file that is in the same directory
    or subdirectory as the file this function is called from.

    Parameters
    ----------
    filename: str
        The name of the file

    subdirs: list of strs, optional
        A list of strings containing any subdirectory names.
        Default: None

    Returns
    -------
    str
        The whole path to the file

    """

    if subdirs is None:
        path = os.path.join(os.path.dirname(__file__), filename)
    elif isinstance(subdirs, list):
        path = os.path.join(os.path.dirname(__file__), *subdirs, filename)
    else:
        msg = ("subdirs must have type list. "
               "If you want a single subdirectory, use subdirs=['data']")
        raise ValueError(msg)

    return path


def calc_mean_from_pdf(x, pdf, dx=None):
    """
    Calculates the mean of a probability density function

    Parameters
    ----------
    x:

    pdf:

    dx:

    """
    if dx is None:
        # If no dx is provided assume they are linearly spaced
        dx = (x[-1] - x[0]) / len(x)

    return np.sum(pdf * x * dx)

def calc_variance_from_pdf(x, pdf, dx=None):
    """
    Calculates the variance from a probability density
    function.

    Parameters
    ----------
    x:

    pdf:

    dx: optional
        Default: None
    """

    if dx is None:
        # If no dx is provided assume they are linearly spaced
        dx = (x[-1] - x[0]) / len(x)

    mean = calc_mean_from_pdf(x, pdf, dx)

    return np.sum(pdf * dx * (x - mean)**2)

def calc_std_from_pdf(x, pdf, dx=None):
    """
    Calculates the standard deviation from a probability
    density function.

    Parameters
    ----------
    x

    pdf

    dx: optional
        Default: None
    """
    if dx is None:
        # If no dx is provided assume they are linearly spaced
        dx = (x[-1] - x[0]) / len(x)

    return np.sqrt(calc_variance_from_pdf(x, pdf, dx))


def calc_z_from_pdf_percentile(x, pdf, percentile):
    """
    """
    cumsum = np.cumsum(pdf)
    normed_cumsum = cumsum / cumsum[-1]

    interpolated_cumsum = interpolate.interp1d(normed_cumsum, x)



    return(interpolated_cumsum(percentile))




def calc_median_from_pdf(x, pdf):
    """
    Calc
    """

    return calc_z_from_pdf_percentile(x, pdf, percentile=0.5)




def normalise_to_pdf(hist, bin_widths):
    """
    """
    if np.sum(hist) < 1e-16:
        pdf = np.zeros(len(hist))
    else:
        pdf = hist/bin_widths/np.sum(hist)

    return pdf


def linear_interpolate_pdfs(sample, xvals, pdfs):
    """
    """
    x1, x2 = xvals
    pdf1, pdf2 = pdfs

    grad = (pdf2 - pdf1) / (x2 - x1)
    dist = sample - x1

    return grad * dist + pdf1



def sigma_to_pdf_percentiles(sigma):
    """


    Parameters
    ----------
    sigma: [1, 2, 3, 4, 5]


    Returns
    -------
    float
        Lower
    float
        Higher
    """

    std = int(sigma)
    std_prop = {
        1: 0.682689492,
        2: 0.954499736,
        3: 0.997300204,
        4: 0.99993666,
        5: 0.999999426697,
    }

    std_limits = {
        1: ((1 - std_prop[1]) / 2, (1 + std_prop[1]) / 2),
        2: ((1 - std_prop[2]) / 2, (1 + std_prop[2]) / 2),
        3: ((1 - std_prop[3]) / 2, (1 + std_prop[3]) / 2),
        4: ((1 - std_prop[4]) / 2, (1 + std_prop[4]) / 2),
        5: ((1 - std_prop[5]) / 2, (1 + std_prop[5]) / 2),
    }

    return std_limits[std]

def redshift_prior(zbins, prior="uniform"):
    """
    """


    available_priors = [
        "uniform",
        "volume",
    ]

    if prior not in available_priors:
        msg = ("'{}' is not in the list of available priors".format(prior))
        raise ValueError(msg)


    if prior == "uniform":
        Pz = np.ones_like(zbins)

    elif prior == "volume":
        msg = "The volume dependent prior has not been implimented yet"
        raise NotImplementedError(msg)


    return Pz