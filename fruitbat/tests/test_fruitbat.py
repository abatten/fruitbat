from fruitbat import Frb
from astropy.coordinates import SkyCoord
from astropy import units as u
import pytest


class Test_Frb_Class(object):

    # Create FRB objects for testing
    frb = Frb("Utmost1", dm=1000, dm_uncert=0.0)
    frb_raj_decj = Frb('Utmost1', dm=1000, raj="11:05:50.0", decj="-8:34:12.0")
    frb_gl_gb = Frb('Utmost1', dm=1000, gl="30.5", gb="-60.2") 

    # Test that methods returns the correct value for DM=1000 and planck2018
    def test_methods(self):
        methods = {
                   "ioka2003" : 0.808948679473659,
                   "inoue2004": 0.9838829149238645,
                   "zhang2018": 1.1092684927893448,
                }
    

        for method in methods.keys():
            z = self.frb.calc_redshift(method=method, cosmology="planck2018")
            assert z == methods[method], "Test failed for: {}".format(method)
   

    def test_frb_calc_skycoords_raj_decj(self):
        ra_str = "11:05:50.0"
        dec_str = "-8:34:12.0"
    
        skycoords = self.frb_raj_decj.calc_skycoords()
        test_skycoords = SkyCoord(ra_str, dec_str, frame="icrs", 
                                   unit=(u.hourangle, u.deg))

        ra, dec = skycoords.ra.value, skycoords.dec.value
        test_ra, test_dec = test_skycoords.ra.value, test_skycoords.dec.value
        assert (ra, dec) == (test_ra, test_dec)
    
    
    def test_frb_calc_skycoords_gl_gb(self):
        gl_str = "30.5"
        gb_str = "-60.2"
    
        skycoords = self.frb_gl_gb.calc_skycoords()
        test_skycoords = SkyCoord(gl_str, gb_str, frame="galactic", unit=u.deg)

        gl, gb = skycoords.l.value, skycoords.b.value
        test_gl, test_gb = test_skycoords.l.value, test_skycoords.b.value
        assert (gl, gb) == (test_gl, test_gb)
    
    
    def test_invalid_method(self):
        invalid_method = "jacqui1992" 
        with pytest.raises(ValueError):
            self.frb.calc_redshift(invalid_method)
            
    
    def test_invalid_cosmology(self):
        invalid_cosmology = "cosmos_1964"
        with pytest.raises(ValueError):
            self.frb.calc_redshift(invalid_cosmology)
