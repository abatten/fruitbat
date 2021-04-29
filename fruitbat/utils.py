"""
The collection of utility functions for Fruitbat.
"""
from __future__ import print_function, absolute_import, division

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
