"""
test_fruitbat_regression.py

This test suite if for testing previously known bugs.
"""
import numpy as np

from fruitbat import Frb, methods, table



def test_create_custom_method():
    """
    This bug was found when creating a custon method,
    using calc_redshift wouldn't return any value.
    """

    def simple_dm(z):
        """ A simple DM-z scaling relation"""
        dm = 1200 * z
        return dm

    methods.add_method("pytest_output_simple_dm", simple_dm)
    table.create("pytest_output_simple_dm")
    f = Frb(1200)
    print(f.calc_redshift(method="pytest_output_simple_dm"))
    assert np.isclose(f.calc_redshift(method="pytest_output_simple_dm"), 1.0)
