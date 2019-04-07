"""
Utility functions for Fruitbat
"""
from __future__ import print_function, absolute_import, division
import os
from six import PY3, PY2

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import astropy.constants as const
import astropy.units as u

__all__ = ["create_lookup_table", "load_lookup_table"]

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
