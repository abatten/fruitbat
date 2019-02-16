from fruitbat import Frb


def test_frb_get_redshift_batten2019():
    f = Frb('Utmost1', 300, 1)
    assert f.calc_redshift(method='batten2019') == 12.0


def test_frb_get_redshift_ioka2003():
    f = Frb('Utmost1', 300, 1)
    assert f.calc_redshift(method='ioka2003') == 100.0


def test_frb_get_redshift_inoue2004():
    f = Frb('Utmost1', 1000, 0)
    assert f.calc_redshift(method='inoue2004', 
        cosmology='planck2018+bao') == (1.1961043734381342)
