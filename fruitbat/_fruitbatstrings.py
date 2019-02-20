"""
Contains a collection of strings that appear many times in the fruitbat
documentation.
"""


def _docstr_sub(*args, **kwargs):
    if args and kwargs:
        raise InputError("Only positional or keyword arguments are allowed")
    else:
        params = args or kwargs

    def do_sub(docstr):
        if docstr.__doc__:
            docstr.__doc__ = docstr.__doc__ % params
        else:
            raise InputError("Function has no docstring for substitution")

        return(docstr)
    return(do_sub)


_methods_doc = "``'batten2019'``, ``'zhang2018'``, ``'inoue2004'``, ``'ioka2003'``"

_cosmo_doc = "``'wmap2013'``, ``planck2013``, ``'planck2015'``, ``planck2018``, ``eagle``"

_dm_units_doc = ":math:`\\rm{pc\\ cm^{-3}}`"
