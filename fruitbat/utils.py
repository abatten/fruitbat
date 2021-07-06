"""
The collection of utility functions for Fruitbat.
"""
import os
import numpy as np
from scipy import interpolate

__all__ = ["check_keys_in_dict", "check_type", "calc_mean_from_pdf", 
    "calc_median_from_pdf", "calc_std_from_pdf", "calc_variance_from_pdf", 
    "calc_z_from_pdf_percentile", "normalise_to_pdf", "redshift_prior", 
    "sigma_to_pdf_percentiles", ]


def check_type(value_name, value, dtype, desire=True):
    """
    Checks the type of a variable and raises an error if not the desired type.

    Parameters
    ----------
    value_name : str
        The name of the variable that will be printed in the error message.

    value :
        The value of the variable

    dtype : dtype
        The data type to compare with isinstance

    desire : boolean, optional
        If `desire = True`, then the error will be raised if value does not
        have a data type of `dtype`. If `desire = False`, then the error will
        be raised if value does have a data type of `dtype`.

    Returns
    -------
    None

    """

    if isinstance(value, dtype) is not desire:

        # Change the error message depending on if we did or did
        # not want a specific data type
        if desire:
            msg_add_in = "have"
        elif not desire:
            msg_add_in = "not have"

        msg = ("The value of {0} should {3} type: {1}. "
               "Instead type({0}) = {2}".format(value_name, dtype, type(value), msg_add_in))

        raise ValueError(msg)

    else:
        pass

def check_keys_in_dict(dictionary, keys):
    """
    Checks that a list of keys exist in a dictionary.

    Parameters
    ----------
    dictionary : dict
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
    filename : str
        The name of the file

    subdirs : list of strs, optional
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
    x : np.ndarray
        The x values.

    pdf : np.ndarray
        The value of the PDF at x.

    dx : np.ndarray or None, optional
        The spacing between the x bins. 
        If `None`, then the bins are assumed to be linearly spaced.

    Returns
    -------
    mean : float
        The mean of the PDF.

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
    x : np.ndarray
        The x values.

    pdf : np.ndarray
        The value of the PDF at x.

    dx : np.ndarray or None, optional
        The spacing between the x bins. 
        If `None`, then the bins are assumed to be linearly spaced.

    Returns
    -------
    variance : float
        The variance of the PDF.

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
    x : np.ndarray
        The x values.

    pdf : np.ndarray
        The value of the PDF at x.

    dx : np.ndarray or None, optional
        The spacing between the x bins. 
        If `None`, then the bins are assumed to be linearly spaced.

    Returns
    -------
    std : float
        The standard deviation of the PDF.

    """
    if dx is None:
        # If no dx is provided assume they are linearly spaced
        dx = (x[-1] - x[0]) / len(x)

    return np.sqrt(calc_variance_from_pdf(x, pdf, dx))


def calc_z_from_pdf_percentile(x, pdf, percentile):
    """


    Parameters
    ----------
    x : np.ndarray
        The x values of the PDF.

    pdf : np.ndarray
        The value of the PDF at x.

    percentile : float
        The percentile of the PDF.

    Returns
    -------
    redshift : float
        The redshift at the given percentile.

    """
    cumsum = np.cumsum(pdf)
    normed_cumsum = cumsum / cumsum[-1]
    interpolated_cumsum = interpolate.interp1d(normed_cumsum, x)
    return interpolated_cumsum(percentile)




def calc_median_from_pdf(x, pdf):
    """
    Calculates the median of a PDF.

    Parameters
    ----------
    x : np.ndarray
        The x values.

    pdf : np.ndarray
        The value of the PDF at x.

    Returns
    -------
    median: float
        The median of the PDF.

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

    Parameters
    ----------
    sample

    xvals: 

    Returns
    -------
    PDF: np.ndarray
        The PDF at sample.
    """
    x1, x2 = xvals
    pdf1, pdf2 = pdfs

    grad = (pdf2 - pdf1) / (x2 - x1)
    dist = sample - x1

    return grad * dist + pdf1



def sigma_to_pdf_percentiles(sigma):
    """
    Looks up the percentile range of Gaussian for a given
    standard deviation.

    Parameters
    ----------
    sigma: [1, 2, 3, 4, 5]
        The standard deviation to calculate a percentile.

    Returns
    -------
    Lower: float
        The lower percentile
    Higher: float
        The higher percentile

    Example
    -------
    >>> sigma_to_pdf_percentiles(1)
    (0.158655254, 0.841344746)

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
