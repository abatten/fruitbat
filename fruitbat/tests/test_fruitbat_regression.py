import os
from glob import glob

import numpy as np

import pytest
import pytest_mpl

from astropy.coordinates import SkyCoord
from astropy import units as u
import pyymw16 as ymw16

from fruitbat import Frb, utils, cosmologies, methods, table, plot, catalogue


# THis bug was found when creating a custon method,
# using calc_redshift wouldn't return any value.
def test_create_custom_method():

    def simple_dm(z):
        dm = 1200 * z
        return dm

    methods.add_method("pytest_output_simple_dm", simple_dm)
    table.create("pytest_output_simple_dm")
    f = Frb(1200)
    print(f.calc_redshift(method="pytest_output_simple_dm"))
    assert np.isclose(f.calc_redshift(method="pytest_output_simple_dm"), 1.0)