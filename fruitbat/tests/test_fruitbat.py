import os

import numpy as np
from glob import glob
import pytest

from astropy.coordinates import SkyCoord
from astropy import units as u
import pyymw16 as ymw16

import fruitbat
from fruitbat import Frb
from fruitbat import utils
from fruitbat import cosmology


class TestFrbClass:

    # Create FRB objects for testing
    frb = Frb(dm=1000, dm_excess=1000)
    frb_raj_decj = Frb(dm=1000, raj="11:05:50.0", decj="-8:34:12.0")
    frb_gl_gb = Frb(dm=1000, gl="30.5", gb="-60.2")
    frb_w_s = Frb(dm=1000, w_obs=30.0, s_peak_obs=20.0)
    frb_host_known = Frb(dm=1000, dm_excess=900, z_host=1.0, dm_host_loc=200)
    frb_dm_host_0 = Frb(dm=1000, dm_excess=900, z_host=1.0)
    frb_dm_host_est = Frb(dm=1100, dm_host_est=100)

    # Test that methods returns the correct value for DM=1000 and planck2018
    def test_methods(self):
        methods = {
            "ioka2003": 0.8089,
            "inoue2004": 0.9838,
            "zhang2018": 1.1092,
        }

        for method in methods.keys():
            z = self.frb.calc_redshift(method=method, cosmology="Planck18")
            assert np.isclose(z, methods[method], rtol=1e-3), "Failed: {}".format(method)

    # Test that a ValueError is raised when an invalid method is given.
    def test_invalid_method(self):
        invalid_method = "jacqui1992"
        with pytest.raises(ValueError):
            self.frb.calc_redshift(method=invalid_method, cosmology="Planck18")

    # Test that a ValueError is raised when an invalid cosmology is given.
    def test_invalid_cosmology(self):
        invalid_cosmology = "cosmos_1964"
        with pytest.raises(ValueError):
            self.frb.calc_redshift(cosmology=invalid_cosmology)

    # Test raises error on dispersion measure less than zero
    def test_frb_negative_dm(self):
        with pytest.raises(ValueError):
            neg_dm = Frb(dm=-1000)

    # Test that the skycoords are calculated correctly when given raj and decj
    def test_frb_calc_skycoords_raj_decj(self):
        ra_str = "11:05:50.0"
        dec_str = "-8:34:12.0"

        skycoords = self.frb_raj_decj.calc_skycoords()
        test_skycoords = SkyCoord(ra_str, dec_str, frame="icrs",
                                  unit=(u.hourangle, u.deg))

        ra, dec = skycoords.ra.value, skycoords.dec.value
        test_ra, test_dec = test_skycoords.ra.value, test_skycoords.dec.value
        assert np.isclose((ra, dec), (test_ra, test_dec)).all()

    # Test that the skycoords are calculated correctly when given gl and gb
    def test_frb_calc_skycoords_gl_gb(self):
        gl_str = "30.5"
        gb_str = "-60.2"

        skycoords = self.frb_gl_gb.calc_skycoords()
        test_skycoords = SkyCoord(gl_str, gb_str, frame="galactic", unit=u.deg)

        gl, gb = skycoords.galactic.l.value, skycoords.galactic.b.value
        test_gl, test_gb = test_skycoords.l.value, test_skycoords.b.value
        assert np.isclose((gl, gb), (test_gl, test_gb)).all()

    # Test that calc_skycoords raises an error if no coords are given
    def test_frb_calc_skycoords_no_coords(self):
        with pytest.raises(ValueError):
            self.frb.calc_skycoords()

    # Test f_obs is calculated correctly when given w_obs and s_peak_obs.
    def test_frb_calc_f_obs(self):
        f_obs = self.frb_w_s.calc_f_obs()
        assert np.isclose(f_obs, 600.0)

    # Test calc_f_obs raises a ValueError if w_obs and s_peak_obs are None.
    def test_frb_calc_f_obs_raise_error(self):
        with pytest.raises(ValueError):
            self.frb.calc_f_obs()

    # Test calc_dm_igm calculates the dm_igm correctly for a known host.
    def test_frb_calc_dm_igm(self):
        dm_igm = self.frb_host_known.calc_dm_igm()
        assert np.isclose(dm_igm, 800.0)

    # Test calc_dm_igm raises ValueError when z is None.
    def test_frb_calc_dm_igm_z_none(self):
        with pytest.raises(ValueError):
            self.frb_w_s.calc_dm_igm()

    # Test calc_dm_igm raises ValueError when dm_host is 0.0 and z is not None.
    def test_frb_calc_dm_igm_dm_host_zero(self):
        with pytest.raises(ValueError):
            self.frb_dm_host_0.calc_dm_igm()

    # Test calc_redshift with subract_host
    def test_frb_calc_redshift_subtract_host(self):
        dm_1 = self.frb_dm_host_est.calc_redshift(subtract_host=True)
        dm_2 = self.frb.calc_redshift()
        assert np.isclose(dm_1, dm_2)
        
    # Test that calc_redshift will raise error if subtract_host is not a bool
    def test_frb_subtract_host_not_bool(self):
        with pytest.raises(ValueError):
            self.frb_dm_host_est.calc_redshift(subtract_host="yes") 

    # Test calc_dm_galaxy calculates dm_galaxy correctly for given coordinates.
    def test_frb_calc_dm_galaxy(self):
        dm_galaxy = self.frb_raj_decj.calc_dm_galaxy()
        dm_pymw16, t_sc_pymw16 = ymw16.dist_to_dm(
            self.frb_raj_decj.skycoords.galactic.l, 
            self.frb_raj_decj.skycoords.galactic.b, 25000)
        assert np.isclose(dm_galaxy, dm_pymw16.value, atol=1e-4)

    # Test calc_dm_galaxy raises a ValueError when no coordinates are given
    def test_frb_cal_dm_galaxy_no_coords(self):
        with pytest.raises(ValueError):
            self.frb.calc_dm_galaxy(model="ymw16")
        
    # Test that the FRB __repr__ is printed
    def test_frb__repr__(self):
        print(self.frb)

    # Test all methods and properties get values and print correctly
    def test_frb_attrs(self):
        for d in dir(self.frb):
            attr = getattr(self.frb, d)
            print(attr)

def test_create_tables():
    method_list = ["ioka2003", "inoue2004", "zhang2018"]
    cosmology_list = fruitbat.cosmology.builtin()

    for method in method_list:
        for key in cosmology_list:
            cosmo = fruitbat.cosmology.builtin()[key]
            outfile_name = "_".join(["pytest_output", method, key])
            fruitbat.utils.create_lookup_table(outfile_name, method=method,
                                               cosmology=cosmo, zmin=0, zmax=3,
                                               num_samples=1000)

    # Remove the files at end of test
    test_files = glob("pytest_output_*")
    for file in test_files:
        os.remove(file)


def test_create_cosmology():

    # Test FlatLambdaCDM
    FlatLambdaCDM_params = {'H0': 67, 'Om0': 0.3, 'flat': True}
    FlatLambdaCDM = cosmology.create_cosmology(FlatLambdaCDM_params)

    # Test FlatwCDM
    FlatwCDM_params = {'H0': 67, 'Om0': 0.3, 'flat': True, 'w0': 0.9}
    FlatwCDM = cosmology.create_cosmology(FlatwCDM_params)

    # Test LambdaCDM
    LambdaCDM_params = {'H0': 67, 'Om0': 0.3, 'Ode0': 0.8, 'flat': False}
    LambdaCDM = cosmology.create_cosmology(LambdaCDM_params)

    # Test wCDM
    wCDM_params = {'H0': 67, 'Om0': 0.3, 'Ode0': 0.8, 'flat': False, 'w0': 0.9}
    wCDM = cosmology.create_cosmology(wCDM_params)



class Test_fz_integrand:

    # Create default cosmology
    cosmo = fruitbat.cosmology.create_cosmology()
    cosmo_w0 = cosmology.create_cosmology({'w0': 1})

    # Test _fz_integrand correctly computes for z = 0
    def test_fz_integrand_z0(self):
        fz = utils._fz_integrand(0, self.cosmo)
        assert np.isclose(fz, 1.0)

    # Test _fz_integrand correctly computes for z = 2
    def test_fz_integrand_z2(self):
        fz = utils._fz_integrand(2, self.cosmo)
        assert np.isclose(fz, 1.011299)

    def test_fz_integrand_w09_z1(self):
        fz = utils._fz_integrand(1, self.cosmo_w0)
        assert np.isclose(fz, 0.291111)


# Test _check_keys_in_dict raises a KeyError when dict is missing keys
def test_check_keys_in_dict_missing():
    required_keys = ["key1", "key2"]
    dictionary = {"key1": 1, "otherkey": 2}
    with pytest.raises(KeyError):
        utils._check_keys_in_dict(dictionary, required_keys)

def test_check_keys_in_dict_all():
    required_keys = ["key1", "key2"]
    dictionary = {"key1": 1, "key2": 2}
    result = utils._check_keys_in_dict(dictionary, required_keys)
    assert result == True
