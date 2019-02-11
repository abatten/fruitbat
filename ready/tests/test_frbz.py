from ready import Frb


def test_frb_get_redshift_batten2019():
    f = Frb('Utmost1', 300, 1)
    assert f.calc_redshift() == (12.0, 4.0)


def test_frb_get_redshift_ioka2003():
    f = Frb('Utmost1', 300, 1)
    assert f.calc_redshift(method='ioka2003') == (100.0, 5.0)


def test_frb_get_redshift_inoue2004():
    f = Frb('Utmost1', 300, 1)
    assert f.calc_redshift(method='inoue2004') == (50.0, 2.0)
