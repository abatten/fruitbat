from fruitbat import Frb
from astropy.coordinates import SkyCoord
from astropy import units as u


def test_methods():
    """
    Test that the different methods calculate the correct redshift for
    planck2018 cosmology with a DM of 1000 pc cm**-3 
    """

    methods = {"batten2019": 12.0,
               "ioka2003" : 100.0,
               "inoue2004": 1.1997438393991349
            }

    frb = Frb("Utmost1", dm=1000, dm_uncert=0.0)
    for method in methods.keys():
        z = frb.calc_redshift(method=method, cosmology="planck2018")
        assert z == methods[method], "Test failed for: {}".format(method)


#def test_frb_get_redshift_batten2019():
#    f = Frb('Utmost1', 300, 1)
#    assert f.calc_redshift(method='batten2019') == 12.0


#def test_frb_get_redshift_ioka2003():
#    f = Frb('Utmost1', 300, 1)
#    assert f.calc_redshift(method='ioka2003') == 100.0


#def test_frb_get_redshift_inoue2004():
#    f = Frb('Utmost1', 1000, 0)
#    assert f.calc_redshift(method='inoue2004', 
#        cosmology='planck2018+bao') == (1.1961043734381342)


def test_frb_calc_skycoords_raj_decj():
    ra_str = "11:05:50.0"
    dec_str = "-12:34:12.0"

    f = Frb('Utmost1', dm=1000, raj=ra_str, decj=dec_str)
    skycoords = f.calc_skycoords()
    test_skycoords = SkyCoord(ra_str, dec_str, frame="icrs", 
                               unit=(u.hourangle, u.deg))
    ra, dec = skycoords.ra.value, skycoords.dec.value
    test_ra, test_dec = test_skycoords.ra.value, test_skycoords.dec.value
    assert (ra, dec) == (test_ra, test_dec)


def test_frb_calc_skycoords_gl_gb():
    gl_str = "30.5"
    gb_str = "-60.2"

    f = Frb('Utmost1', dm=1000, gl=gl_str, gb=gb_str)
    skycoords = f.calc_skycoords()
    test_skycoords = SkyCoord(gl_str, gb_str, frame="galactic", unit=u.deg)
    gl, gb = skycoords.l.value, skycoords.b.value
    test_gl, test_gb = test_skycoords.l.value, test_skycoords.b.value
    assert (gl, gb) == (test_gl, test_gb)

