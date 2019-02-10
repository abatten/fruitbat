from frbz import redshift_estimation as re
from frbz import Frb

def test_dm_to_redshift_no_err():
    z = re.dm_to_redshift(9.0)
    assert z == (2.0, 0.0)

def test_dm_to_redshift_with_err():
    z = re.dm_to_redshift(9.0, 4.0)
    assert z == (2.0, 1.0)


def test_frb_get_redshift_batten2019():
    f = Frb.Frb('Utmost1', 300, 1)
    assert f.calc_redshift() == (12.0, 4.0)
