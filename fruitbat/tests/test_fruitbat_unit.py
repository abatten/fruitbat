import os
from glob import glob
from textwrap import dedent

import numpy as np

import pytest
import pytest_mpl

import astropy
from astropy.coordinates import SkyCoord
from astropy import units as u
import pygedm

import fruitbat
from fruitbat import Frb, utils, cosmologies, methods, table, plot, catalogue


def test_get_bibtex_is_correct():

    ads_bibtex = dedent(
    r"""
        @ARTICLE{2019JOSS....4.1399B,
               author = {{Batten}, Adam},
                title = "{Fruitbat: A Python Package for Estimating Redshifts of Fast Radio Bursts}",
              journal = {The Journal of Open Source Software},
             keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - High Energy Astrophysical Phenomena},
                 year = "2019",
                month = "May",
               volume = {4},
               number = {37},
                pages = {1399},
                  doi = {10.21105/joss.01399},
        archivePrefix = {arXiv},
               eprint = {1905.04294},
         primaryClass = {astro-ph.IM},
               adsurl = {https://ui.adsabs.harvard.edu/abs/2019JOSS....4.1399B},
              adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }
    """
    )
    assert fruitbat.get_bibtex() == ads_bibtex.strip()

def test_cite_is_get_bibtex():
    assert fruitbat.__cite__() == fruitbat.get_bibtex()


class TestFrbClass:

    # Create FRB objects for testing
    frb = Frb(dm=1000, dm_excess=1000, name='simple_frb')
    frb_raj_decj = Frb(dm=1000, raj="11:05:50.0", decj="-8:34:12.0")
    frb_gl_gb = Frb(dm=1000, gl="30.5", gb="-60.2")
    frb_w_s = Frb(dm=1000, width=30.0, peak_flux=20.0)
    frb_host_known = Frb(dm=1000, dm_excess=900, z_host=1.0, dm_host_loc=200)
    frb_dm_host_0 = Frb(dm=1000, dm_excess=900, z_host=1.0)
    frb_dm_host_est = Frb(dm=1100, dm_host_est=100)
    frb_energy = Frb(dm=1000, obs_bandwidth=400, width=1, peak_flux=2)
    frb_energy_freq = Frb(dm=1000, obs_freq_central=0.4, width=1, peak_flux=2)
    frb_utc = Frb(dm=1000, utc="1999-01-01T00:00:00.000")
    frb_with_units = Frb(dm=1000, obs_bandwidth=400*u.MHz)
    frb_fluence = Frb(dm=1000, fluence=2)


    def test_frb_dm_assignment_int(self):
        test_frb = Frb(1000)
        assert test_frb.dm == 1000.0 * u.pc * u.cm**-3

    def test_frb_dm_assignment_float(self):
        test_frb = Frb(1000.0)
        assert test_frb.dm == 1000.0 * u.pc * u.cm**-3

    def test_frb_dm_assignment_string(self):
        with pytest.raises(ValueError):
            test_frb = Frb("1000")

    def test_frb_gl_gb_assignment_float(self):
        test_frb = Frb(1000, gl=10, gb=20)
        assert (test_frb.gl.value, test_frb.gb.value) == (10, 20)

    def test_frb_gl_gb_assignment_string(self):
        test_frb = Frb(1000, gl="10", gb="20")
        assert (test_frb.gl.value, test_frb.gb.value) == (10, 20)


    # Test that methods returns the correct value for DM=1000 and planck2018
    def test_methods(self):
        methods = {
            "Ioka2003": 0.80856155,
            "Inoue2004": 0.98344417,
            "Zhang2018": 1.10879646
        }

        for method in methods.keys():
            z = self.frb.calc_redshift(method=method, cosmology="Planck18")
            assert np.isclose(z.value, methods[method]), "Fail: {}".format(method)

    # Test that a ValueError is raised when an invalid method is given.
    def test_invalid_method(self):
        invalid_method = "jacqui1992"
        with pytest.raises(ValueError):
            self.frb.calc_redshift(method=invalid_method, cosmology="Planck18")

    # Test that a ValueError is raised when an invalid cosmology is given.
    def test_invalid_cosmology(self):
        invalid_cosmology = "cosmos_1964"
        with pytest.raises(ValueError):
            self.frb.calc_redshift(method="Ioka2003", cosmology=invalid_cosmology)

    # Test raises error on dispersion measure less than zero
    def test_frb_negative_dm(self):
        with pytest.raises(ValueError):
            Frb(dm=-1000)

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

    # Test fluence is calculated correctly when given width and peak_flux.
    def test_frb_calc_fluence(self):
        fluence = self.frb_w_s.calc_fluence()
        assert np.isclose(fluence.value, 600.0)

    # Test calc_fluence raises a ValueError if width and peak_flux are None.
    def test_frb_calc_fluence_raise_error(self):
        with pytest.raises(ValueError):
            self.frb.calc_fluence()

    # Test calc_dm_igm calculates the dm_igm correctly for a known host.
    def test_frb_calc_dm_igm(self):
        dm_igm = self.frb_host_known.calc_dm_igm()
        assert np.isclose(dm_igm.value, 800.0)

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
        dm_1 = self.frb_dm_host_est.calc_redshift(method="Ioka2003", subtract_host=True)
        dm_2 = self.frb.calc_redshift(method="Ioka2003")
        assert np.isclose(dm_1, dm_2)

    # Test that calc_redshift will raise error if subtract_host is not a bool
    def test_frb_subtract_host_not_bool(self):
        with pytest.raises(ValueError):
            self.frb_dm_host_est.calc_redshift(subtract_host="yes")

    # Test calc_dm_galaxy calculates dm_galaxy correctly for given coordinates.
    def test_frb_calc_dm_galaxy_ymw16(self):
        dm_galaxy = self.frb_raj_decj.calc_dm_galaxy("YMW16")
        dm_pymw16, t_sc_pymw16 = pygedm.dist_to_dm(
            self.frb_raj_decj.skycoords.galactic.l,
            self.frb_raj_decj.skycoords.galactic.b, 25000, method="YMW16")
        assert np.isclose(dm_galaxy.value, dm_pymw16.value)


    # Test calc_dm_galaxy calculates dm_galaxy correctly for given coordinates.
#    def test_frb_calc_dm_galaxy_ne2001(self):
#        dm_galaxy = self.frb_raj_decj.calc_dm_galaxy("NE2001")
#        dm_ne2001, t_sc_ne2001 = pygedm.dist_to_dm(
#            self.frb_raj_decj.skycoords.galactic.l,
#            self.frb_raj_decj.skycoords.galactic.b, 25000, method="NE2001")
#        assert np.isclose(dm_galaxy.value, dm_ne2001.value)

    # Test calc_dm_galaxy raises a ValueError when no coordinates are given
    def test_frb_cal_dm_galaxy_no_coords(self):
        with pytest.raises(ValueError):
            self.frb.calc_dm_galaxy(model="ymw16")

    def test_frb_calc_lum_dist_without_z(self):
        with pytest.raises(ValueError):
            self.frb.z = None
            self.frb.calc_luminosity_distance()

    # Test calc_energy calculates the energy of an FRB
    def test_frb_calc_energy_bandwidth(self):
        self.frb_energy.calc_redshift(method="Inoue2004")
        energy = self.frb_energy.calc_energy(use_bandwidth=True)
        assert np.isclose(energy.value, 2.13256754066293e+40)

    def test_frb_calc_energy_frequency(self):
        self.frb_energy_freq.calc_redshift(method="Inoue2004")
        energy = self.frb_energy_freq.calc_energy()
        assert np.isclose(energy.value, 2.13256754066293e+37)

    def test_frb_calc_energy_no_fluence(self):
        with pytest.raises(ValueError):
            self.frb.calc_redshift(method="Inoue2004")
            self.frb.calc_energy(use_bandwidth=True)

    def test_frb_calc_energy_no_bandwidth(self):
        with pytest.raises(ValueError):
            self.frb_fluence.calc_redshift(method="Inoue2004")
            self.frb_fluence.calc_energy(use_bandwidth=True)

    def test_frb_calc_energy_no_frequency(self):
        with pytest.raises(ValueError):
            self.frb_energy.calc_redshift(method="Inoue2004")
            self.frb_energy.calc_energy()

    def test_frb_calc_luminosity_bandwidth(self):
        self.frb_energy.calc_redshift(method="Inoue2004")
        lum = self.frb_energy.calc_luminosity(use_bandwidth=True)
        assert np.isclose(lum.value, 4.229828665e+43)

    def test_frb_calc_luminosity_frequency(self):
        self.frb_energy_freq.calc_redshift(method="Inoue2004")
        lum = self.frb_energy_freq.calc_luminosity()
        assert np.isclose(lum.value, 4.2298286655e+40)

    def test_frb_calc_luminosity_no_frequency(self):
        with pytest.raises(ValueError):
            self.frb_energy.calc_redshift(method="Inoue2004")
            self.frb_energy.calc_luminosity()

    def test_frb_calc_comoving_distance(self):
        self.frb.calc_redshift(method="Inoue2004")
        dist = self.frb.calc_comoving_distance()
        assert np.isclose(dist.value, 3351.51321266)

    def test_frb_pass_wrong_units(self):
        with pytest.raises(ValueError):
            Frb(dm=1000, obs_bandwidth=400*u.m)

    # Test that the FRB __repr__ is printed
    def test_frb__repr__(self):
        print(self.frb)

    # Test all methods and properties get values and print
    def test_frb_attrs(self):
        for d in dir(self.frb):
            attr = getattr(self.frb, d)
            print(attr)


def test_create_FlatLambdaCDM_cosmology():
    # Test FlatLambdaCDM
    FlatLambdaCDM_params = {'H0': 67, 'Om0': 0.3, 'flat': True}
    flcdm = cosmologies.create_cosmology(FlatLambdaCDM_params)
    assert type(flcdm) == astropy.cosmology.core.FlatLambdaCDM

def test_create_FlatwCDM_cosmology():
    # Test FlatwCDM
    FlatwCDM_params = {'H0': 67, 'Om0': 0.3, 'flat': True, 'w0': 0.9}
    fwcdm = cosmologies.create_cosmology(FlatwCDM_params)
    assert type(fwcdm) == astropy.cosmology.core.FlatwCDM

def test_create_LambdaCDM_cosmology():
    # Test LambdaCDM
    LambdaCDM_params = {'H0': 67, 'Om0': 0.3, 'Ode0': 0.8, 'flat': False}
    lcdm = cosmologies.create_cosmology(LambdaCDM_params)
    assert type(lcdm) == astropy.cosmology.core.LambdaCDM

def test_create_wCDM_cosmology():
    # Test wCDM
    wCDM_params = {'H0': 67, 'Om0': 0.3, 'Ode0': 0.8, 'flat': False, 'w0': 0.9}
    wcdm = cosmologies.create_cosmology(wCDM_params)
    assert type(wcdm) == astropy.cosmology.core.wCDM


class Test_fz_integrand:

    # Create default cosmology
    cosmo = cosmologies.create_cosmology()
    cosmo_w0 = cosmologies.create_cosmology({'w0': 1})

    # Test _fz_integrand correctly computes for z = 0
    def test_fz_integrand_z0(self):
        fz = methods._f_integrand(0, self.cosmo)
        assert np.isclose(fz, 1.0)

    # Test _fz_integrand correctly computes for z = 2
    def test_fz_integrand_z2(self):
        fz = methods._f_integrand(2, self.cosmo)
        assert np.isclose(fz, 1.011299)

    def test_fz_integrand_w1_z1(self):
        fz = methods._f_integrand(1, self.cosmo_w0)
        assert np.isclose(fz, 0.291111)


# Test _check_keys_in_dict raises a KeyError when dict is missing keys
def test_check_keys_in_dict_missing():
    required_keys = ["key1", "key2"]
    dictionary = {"key1": 1, "otherkey": 2}
    with pytest.raises(KeyError):
        utils.check_keys_in_dict(dictionary, required_keys)


def test_check_keys_in_dict_all():
    required_keys = ["key1", "key2"]
    dictionary = {"key1": 1, "key2": 2}
    result = utils.check_keys_in_dict(dictionary, required_keys)
    assert result


class TestAddingMethods:

    def new_method(self, z, cosmo):
        return 1200 * z

    def test_add_method(self):
        methods.add_method("new_method", self.new_method)
        assert "new_method" in methods.available_methods()

    def test_reset_methods(self):
        methods.reset_methods()
        assert "new_method" not in methods.available_methods()


class TestCatalogue:

    def test_create_analysis_catalogue(self):
        catalogue.create_analysis_catalogue("pytest_output_analysis_catalogue")
        assert os.path.exists("pytest_output_analysis_catalogue.csv")


    def test_create_method_catalogue(self):
        catalogue.create_methods_catalogue("pytest_output_methods_catalogue")
        assert os.path.exists("pytest_output_methods_catalogue.csv")


class TestCreateTables:

    def test_create_tables_normal(self):
        method_list = methods.builtin_method_functions()
        cosmology_list = cosmologies.builtin_cosmology_functions()

        # Create a lookup table for each analytic method and cosmology
        for method_key in method_list:
            if method_key in methods.methods_hydrodynamic():
                continue
            for cosmo_key in cosmology_list:
                if cosmo_key == "EAGLE":  # Skip EAGLE since it is the same as Planck13
                    continue
                here = os.getcwd()

                cosmo = cosmologies.builtin_cosmology_functions()[cosmo_key]
                filename = "_".join(["pytest_output", method_key, cosmo_key])
                table.create(method=method_key, filename=filename,
                             cosmo=cosmo, output_dir=here, zmin=0,
                             zmax=20, num_samples=10000)

                # Compare new tables to existing tables for 4 dm values
                pre_calc_table_name = "{}.hdf5".format(method_key)

                pre_calc_table = utils.get_path_to_file_from_here(pre_calc_table_name, subdirs=["data"])
                new_calc_table = "pytest_output_{}_{}.hdf5".format(method_key, cosmo_key)

                #pre_calc = table.load(pre_calc_fn)
                #new_calc = table.load(new_calc_fn, data_dir=here)

                test_dm_list = [0, 100, 1000, 2000]

                for dm in test_dm_list:
                    new_z = table.get_z_from_table(dm, new_calc_table, cosmo_key)
                    pre_z = table.get_z_from_table(dm, pre_calc_table, cosmo_key)
                    assert np.isclose(new_z, pre_z, rtol=1e-03)

    def test_create_table_zhang_figm_free_elec(self):
        cosmo = cosmologies.builtin_cosmology_functions()["Planck18"]
        filename = "_".join(["pytest_output", "Zhang2018",
                             "Planck18", "figm_free_elec"])
        here = os.getcwd()

        table.create(method="Zhang2018", filename=filename, cosmo=cosmo,
                     output_dir=here, f_igm=0.5, free_elec=0.4)

    def test_create_table_zhang_figm_error(self):
        cosmo = cosmologies.builtin_cosmology_functions()["Planck18"]

        with pytest.raises(ValueError):
            table.create(method="Zhang2018", cosmo=cosmo, f_igm=-1)

    def test_create_table_zhang_free_elec_error(self):
        cosmo = cosmologies.builtin_cosmology_functions()["Planck18"]
        filename = "_".join(["pytest_output", "Zhang2018",
                             "Planck18", "free_elec_error"])

        with pytest.raises(ValueError):
            table.create(method="Zhang2018", filename=filename, cosmo=cosmo,
                         free_elec=-1)

    def test_create_table_invalid_method(self):
        with pytest.raises(ValueError):
            table.create(method="Webb1995")


class TestPlots:
    # Test that the method plot creates an output file
    def test_method_plot(self):
        with pytest_mpl.plugin.switch_backend('Agg'):
            plot.method_comparison(filename="pytest_output_method")
        cwd = os.getcwd()
        if not os.path.exists(os.path.join(cwd, "pytest_output_method.png")):
            raise OSError

    # Test that the cosmology plot creates and output file
    def test_cosmology_plot(self):
        with pytest_mpl.plugin.switch_backend('Agg'):
            plot.cosmology_comparison(filename="pytest_output_cosmo")
        cwd = os.getcwd()
        if not os.path.exists(os.path.join(cwd, "pytest_output_cosmo.png")):
            raise OSError

    def test_redshift_pdf_plot(self):
        frb = Frb(510, gl=34,gb=15)
        with pytest_mpl.plugin.switch_backend('Agg'):
            plot.redshift_pdf(frb, filename="pytest_output_pdf", usetex=False)
        cwd = os.getcwd()
        if not os.path.exists(os.path.join(cwd, "pytest_output_cosmo.png")):
            raise OSError  

            
def test_cleanup():
    # Remove the files at end of test
    test_files = glob("*pytest_output*")
    test_files += glob("./data/*pytest_output*")

    for file in test_files:
        os.remove(file)
