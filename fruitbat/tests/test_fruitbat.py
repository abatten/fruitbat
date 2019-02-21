import fruitbat
from fruitbat import Frb
from fruitbat import utils

from astropy.coordinates import SkyCoord
from astropy import units as u
import pytest
import os
from glob import glob

class Test_Frb_Class(object):

    # Create FRB objects for testing
    frb = Frb("Utmost1", dm=1000, dm_excess=1000, dm_uncert=0.0)
    frb_raj_decj = Frb('Utmost2', dm=1000, raj="11:05:50.0", decj="-8:34:12.0")
    frb_gl_gb = Frb('Utmost3', dm=1000, gl="30.5", gb="-60.2") 
    frb_w_s = Frb('Utmost4', dm=1000, w_obs=30.0, s_peak_obs=20.0)

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
   

    # Test that the skycoords are calculated correctly when given raj and decj
    def test_frb_calc_skycoords_raj_decj(self):
        ra_str = "11:05:50.0"
        dec_str = "-8:34:12.0"
    
        skycoords = self.frb_raj_decj.calc_skycoords()
        test_skycoords = SkyCoord(ra_str, dec_str, frame="icrs", 
                                   unit=(u.hourangle, u.deg))

        ra, dec = skycoords.ra.value, skycoords.dec.value
        test_ra, test_dec = test_skycoords.ra.value, test_skycoords.dec.value
        assert (ra, dec) == (test_ra, test_dec)
    

    # Test that the skycoords are calculated correctly when given gl and gb
    def test_frb_calc_skycoords_gl_gb(self):
        gl_str = "30.5"
        gb_str = "-60.2"
    
        skycoords = self.frb_gl_gb.calc_skycoords()
        test_skycoords = SkyCoord(gl_str, gb_str, frame="galactic", unit=u.deg)

        gl, gb = skycoords.l.value, skycoords.b.value
        test_gl, test_gb = test_skycoords.l.value, test_skycoords.b.value
        assert (gl, gb) == (test_gl, test_gb)
    

    # Test that f_obs is calculated correctly when given w_obs and s_peak_obs.
    def test_frb_calc_f_obs(self):
        f_obs = self.frb_w_s.calc_f_obs()
        assert f_obs == 600.0


    # Test that a ValueError is raised when an invalid method is given.
    def test_invalid_method(self):
        invalid_method = "jacqui1992" 
        with pytest.raises(ValueError):
            self.frb.calc_redshift(method=invalid_method)
            
    # Test that a ValueError is raised when an invalid cosmology is given.
    def test_invalid_cosmology(self):
        invalid_cosmology = "cosmos_1964"
        with pytest.raises(ValueError):
            self.frb.calc_redshift(cosmology=invalid_cosmology)


def test_create_tables():
    method_list = ["ioka2003", "inoue2004", "zhang2018"]
    cosmology_list = fruitbat.cosmology.cosmology_keys()
 
    for method in method_list:
        for cosmology in cosmology_list:
            cosmo_params = fruitbat.cosmology.avaliable_cosmologies()[cosmology]
            outfile_name = "_".join(["pytest_output", method, cosmology])
            fruitbat.utils.create_lookup_table(outfile_name, method=method,
                cosmology=cosmo_params, zmin=0, zmax=3, num_samples=1000)

    # Remove the files at end of test
    test_files = glob("pytest_output_*")
    for file in test_files: 
        os.remove(file)



def test_fz_integrand():
    cosmology = {"Omega_m": 0.3, "Omega_L": 0.7}
    fz = utils._fz_integrand(0, cosmology)
    assert fz == 1.0

