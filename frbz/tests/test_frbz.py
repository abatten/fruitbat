from frbz import redshift_estimation as re

def test_dm_to_redshift_no_err():
    z = re.dm_to_redshift(9.0)
    assert z == (2.0, 0.0)

def test_dm_to_redshift_with_err():
    z = re.dm_to_redshift(9.0, 4.0)
    assert z == (2.0, 1.0)
