"""
Contains a collection of strings that appear many times in the fruitbat
documentation.
"""


def docstr_sub(*args, **kwargs):
    params = args or kwargs

    def do_sub(docstr):
        if docstr.__doc__:
            docstr.__doc__ = docstr.__doc__ % params
        else:
            raise InputError("Function has no docstring for substitution")

        return(docstr)
    return(do_sub)


_methods_doc = "``'batten2019'``, ``'zhang2018'``, ``'inoue2004'``, ``'ioka2003'``"

_cosmo_doc = "``'wmap2012'``, ``'planck2015'``, ``planck2018``, ``planck2018+bao``"

_dm_units_doc = ":math:`\\rm{pc\\ cm^{-3}}`"
