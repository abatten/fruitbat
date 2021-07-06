"""
Frb class and method functions.
"""
import numpy as np
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import h5py

from e13tools import docstring_substitute, docstring_copy
import pygedm

from fruitbat import cosmologies, methods, table, utils, plot

__all__ = ["Frb"]

class Frb(object):
    """
    Create a :class:`~Frb` object using the observered properties of a
    FRB including dispersion measure (DM) and its sky coordinates. This
    class utilises various DM-redshift relations as well as the YMW16
    galactic DM model in the analysis.

    A FRB can be defined simply with a single :attr:`dm` value. This
    should be the observed DM prior to subtracting the Milky Way,
    however if you do not provide any additional information, it is
    asssumed that this value is also the :attr:`dm_excess` (DM -
    Galactic DM). You can get the estimated redshift by using the
    method :meth:`calc_redshift`.

    To get a more accurate distance estimate you can account for the
    contribution due to the Milky Way by supplying :attr:`dm_galaxy` or
    by giving the coordinates of the FRB in ICRS (:attr:`raj`,
    :attr:`decj`) or Galactic (:attr:`gl`, :attr:`gb`) coordinates and
    calling the method :meth:`calc_dm_galaxy` before calling
    :meth:`calc_redshift`.


    Parameters
    ----------
    dm : float
        The observed dispersion measure of the FRB. This is without
        Milky Way or host galaxy subtraction. Units: pc cm**-3

    Keyword Arguments
    -----------------
    name : str or None, optional
        The name of the frb object. Default: *None*

    raj : str or None, optional
        The right ascension in J2000 coordinates of the best estimate
        of the FRB position. By default should be given in the ICRS
        frame. An alternative frame can be specified by also including
        the keyword ``frame``. e.g. ``frame='fk5'``. Default: *None*

    decj : str or None, optional
        The declination in J2000 coordinates of the best estimate of
        the FRB position. By default should be given in the ICRS frame.
        An alternative frame can be specified by also including the
        keyword ``frame``. e.g. ``frame='fk5'``. Default: *None*

    gl : str or None, optional
        The Galactic longitude in degrees of the best estimate of the
        FRB position. Default: *None*

    gb : str or None, optional
        The Galactic latitude in degrees of the best estimate of the
        FRB position. Default: *None*

    Other Parameters
    ----------------
    dm_galaxy : float, optional
        The modelled contribution to the FRB DM by electrons in the
        Milky Way. This value is calculated using
        :meth:`calc_dm_galaxy` Units: pc cm**-3 Default: 0.0

    dm_excess : float or None, optional
        The DM excess of the FRB over the estimated Galactic DM. If
        :attr:`dm_excess` is *None*, then :attr:`dm_excess` is
        calculated automatically with :meth:`calc_dm_excess`.
        Units: pc cm**-3 Default: *None*

    z_host : float or None, optional
        The observed redshift of the localised FRB host galaxy.
        Default: *None*

    dm_host_est : float, optional
        The estimated contribution to the measured FRB DM originating
        from the FRB's host galaxy. This value is the amount of DM the
        host galaxy contributes to the observed DM, *not* the DM of the
        host galaxy. Setting this to a non-zero value and setting
        ``subtract_host=True`` when calling :meth:`calc_redshift`
        accounts for the DM contribution due to the host galaxy.
        Units: pc cm**-3 Default: 0.0

    dm_host_loc : float, optional
        The dispersion measure of a localised FRB host galaxy. This
        value is *not* the contribution to the observed DM, but the DM
        at the host galaxy. The observed DM is :attr:`dm_host_loc` but
        attenuated by a factor of (1 + z). To use this, the redshift of
        the host galaxy must be known. Units: pc cm**-3 Default: 0.0

    dm_index : float or None, optional
        The dispersion measure index of the burst :math:`\\alpha` such
        that :math:`\\rm{DM} \\propto \\nu^{-\\alpha}` Default: *None*

    scatt_index : float or None, optional
        The scattering index (:math:`\\beta`) of the FRB pulse. The
        scattering index describes how the width (:math:`\\rm{W}`) of
        the FRB pulse evolves with frequency :math:`\\nu` such that
        :math:`\\rm{W} \\propto \\nu^{-\\beta}`. Default: *None*

    snr : float or None, optional
        The signal-to-noise ratio of the burst.
        Default: *None*

    width : float or None, optional
        The observed width of the pulse obtained by a pulse fitting
        algorithm. Units: ms Default: *None*

    peak_flux : float or None, optional
        The observed peak flux density of the burst. Units: Jy
        Default: *None*

    fluence : float or None, optional
        The observed fluence of the FRB. If :attr:`fluence` is *None*
        and both :attr:`width` and :attr:`peak_flux` are not *None*
        then :attr:`fluence` is automatically calculated with
        :meth:`calc_fluence`. Units: Jy ms Default: *None*

    obs_bandwidth : float or None, optional
        The observing bandwidth. Units: MHz Default: *None*

    obs_freq_central : float or None, optional
        The central observing frequency, Units: MHz Deault: *None*

    utc : str or None, optional
        The UTC time of the FRB Burst. Format should be of the form
        '1999-01-01T00:00:00.000'. Default: *None*


    Example
    -------
    >>> import fruitbat
    >>> FRB = fruitbat.Frb(879, gl="12:31:40.5", gb="3:41:10.0")
    >>> FRB.calc_dm_galaxy()
    >>> FRB.calc_redshift()

    """

    def __init__(self, dm, name=None, raj=None, decj=None, gl=None, gb=None,
                 dm_galaxy=0.0, dm_excess=None, z_host=None, dm_host_est=0.0,
                 dm_host_loc=0.0, dm_index=None, scatt_index=None, snr=None,
                 width=None, peak_flux=None, fluence=None, obs_bandwidth=None,
                 obs_freq_central=None, utc=None, **kwargs):

        self.dm = dm

        self.name = name
        self.raj = raj
        self.decj = decj
        self.gl = gl
        self.gb = gb

        self.dm_galaxy = dm_galaxy

        # Calculate dm_excess from existing parameters if it is not given.
        if dm_excess is None:
            self.calc_dm_excess()
        else:
            self.dm_excess = dm_excess

        self.z_host = z_host
        self.dm_host_est = dm_host_est
        self.dm_host_loc = dm_host_loc
        self.dm_index = dm_index
        self.scatt_index = scatt_index
        self.snr = snr
        self.width = width
        self.peak_flux = peak_flux
        self.fluence = fluence
        self.obs_bandwidth = obs_bandwidth
        self.obs_freq_central = obs_freq_central
        self.utc = utc

        self.z = None
        self.dm_igm = None
        self.cosmology = None
        self.method = None
        self.dm_galaxy_model = None
        self.z_conf_int_lower = None
        self.z_conf_int_upper = None
        self.tau_sc = None

        # Calculate fluence if peak_flux and width are given
        if ((self.fluence is None) and
                (peak_flux is not None and width is not None)):
            self.calc_fluence()

        # Calculate the skycoords of the FRB from (raj, decj) or (gl, gb)
        if (raj and decj) or (gl and gb):
            if 'frame' in kwargs:
                self.skycoords = self.calc_skycoords(frame=kwargs["frame"])
            else:
                self.skycoords = self.calc_skycoords()

            self.raj = self.skycoords.transform_to('icrs').ra
            self.decj = self.skycoords.transform_to('icrs').dec
            self.gl = self.skycoords.transform_to('galactic').l
            self.gb = self.skycoords.transform_to('galactic').b
        else:
            self.skycoords = None

    def __repr__(self):

        frb_repr = []

        var_dict = vars(self)

        for var in var_dict:
            if var[1:] == "skycoords":
                continue
            var_name = var[1:]

            if isinstance(var_dict[var], u.Quantity):
                var_val = var_dict[var].value
            else:
                var_val = var_dict[var]
            frb_repr.append("{}={}".format(var_name, var_val))

        return "Frb({})".format(", ".join(map(str, frb_repr)))

    @docstring_substitute(meth=methods.available_methods(),
                          cosmo=cosmologies.available_cosmologies())
    def calc_redshift(self, method='Batten2021', cosmology="Planck18",
                      subtract_host=False, lookup_table=None):
        """
        Calculate the redshift of the FRB from its :attr:`dm`,
        :attr:`dm_excess` or :attr:`dm_excess` - :attr:`dm_host_est`.

        Parameters
        ----------
        method : str, optional
            The dispersion meausre -redshift relation to use when
            calculating the redshift. Avaliable methods:  %(meth)s.
            Default: 'Inoue2004'

        cosmology : str, optional
            The cosmology to assume when calculating the redshift.
            Avaliable cosmologies: %(cosmo)s. Default: 'Planck18'

        subtract_host : bool, optional
            Subtract :attr:`dm_host_est` from the :attr:`dm_excess`
            before calculating the redshift. This is is used to account
            for the dispersion measure that arises from the FRB host
            galaxy. Default: False

        lookup_table : str or None, optional
            The path to the lookup table file. If ``lookup_table=None``
            a table will attempted to be loaded from the data directory
            based on the method name. Default: *None*

        Returns
        -------
        float
            The redshift of the FRB.

        Notes
        -----
        The methods_ section in the documentation has a description for
        each of the builtin methods.

        The cosmology_ section in the documentation has a list of the
        cosmological parameters for each cosmology

        .. _cosmology: https://fruitbat.readthedocs.io/en/latest/user_guide/method_and_cosmology.html#cosmology
        .. _methods: https://fruitbat.readthedocs.io/en/latest/user_guide/method_and_cosmology.html#methods

        """
        utils.check_type("subtract_host", subtract_host, bool)

        if subtract_host:
            input_dm = self.dm_excess - self.dm_host_est
            print(input_dm)
        else:
            input_dm = self.dm_excess

        # If the user provides a table use that table for estimation.
        if lookup_table is not None:
            raise NotImplementedError("Not made yet!")

        if method not in methods.available_methods():
            err_msg = ("Method '{}' is not a valid method. "
                       "Valid methods are: {}".format(method,
                        methods.available_methods()))
            raise ValueError(err_msg)

        if method in methods.methods_hydrodynamic():
            zvals, pdf, dz = self.calc_redshift_pdf(method=method)
            z_median = utils.calc_median_from_pdf(zvals, pdf)
            self.z = z_median
            self.method = method
            if method == "Batten2020":
                self.cosmology = "EAGLE"

        elif method in methods.methods_analytic():
            if cosmology in cosmologies.available_cosmologies():
                #if method in methods.builtin_method_functions().keys():
                table_name = table.get_table_path(method)
                self.z = table.get_z_from_table(input_dm, table_name, cosmology)
                self.cosmology = cosmology
                self.method = method


            else:
                err_msg = ("Cosmology '{}' is not a valid cosmology. Valid "
                           "cosmologies are: {}".format(cosmology,
                            cosmologies.available_cosmologies()))
                raise ValueError(err_msg)



        else:
            table_name = utils.get_path_to_file_from_here("{}.hdf5".format(method), subdirs=["data"])
            self.z = table.get_z_from_table(input_dm, table_name)
            self.cosmology = cosmology
            self.method = method

        return self.z

    #@docstring_substitute(meth=methods.available_methods(),
    #                      cosmo=cosmologies.available_cosmologies())
    def calc_redshift_pdf(self, method="Batten2021", cosmology="Planck18", prior="uniform", subtract_host=False,
                          lookup_table=None):
        """
        Calc
        """
        if method == "Batten2021":
            filename = utils.get_path_to_file_from_here("Batten2021.hdf5", subdirs=["data"])


        elif method in methods.methods_analytic():
            filename = "{}.hdf5".format(method)
            filename = utils.get_path_to_file_from_here(filename, subdirs=["data"])
            raise NotImplementedError("Getting PDFs from analytic moddels has not been implemented yet!")

        else:
            raise NotImplementedError("This has not been implemented yet!")


        with h5py.File(filename, "r") as data:

            DMzHist = data["DMz_hist"][:]
            redshift_bin_widths = data["Redshift_Bin_Widths"][:]
            redshifts = data["Redshifts_Bin_Edges"][1:]
            DMBins = data["DM_Bin_Edges"][:]

            max_bin_idx = np.where(self.dm_excess.value <= DMBins)[0][0]
            prev_bin_idx = max_bin_idx - 1

            hist_higher_dm = DMzHist[max_bin_idx]
            hist_lower_dm = DMzHist[prev_bin_idx]

            pdf2 = utils.normalise_to_pdf(hist_higher_dm, bin_widths=redshift_bin_widths)
            pdf1 = utils.normalise_to_pdf(hist_lower_dm, bin_widths=redshift_bin_widths)


            DMlower, DMhigher = DMBins[prev_bin_idx], DMBins[max_bin_idx]
            lin_interp_pdf = utils.linear_interpolate_pdfs(self.dm_excess.value,
                (DMlower, DMhigher), (pdf1, pdf2))

            prior = utils.redshift_prior(redshifts, prior=prior)

            z_pdf = lin_interp_pdf * prior
            z_bins = redshifts
            dz = redshift_bin_widths



        self.z_bins = z_bins
        self.z_pdf = z_pdf
        self.dz = dz

        return z_bins, z_pdf, dz


    @docstring_copy(plot.redshift_pdf)
    def plot_redshift_pdf(self, *args, **kwargs):
        return plot.redshift_pdf(self, *args, **kwargs)


    @docstring_substitute(meth=methods.available_methods(),
                          cosmo=cosmologies.available_cosmologies())
    def calc_redshift_conf_int(self, method="Batten2021", sigma=1, scatter_percentage=0, **calc_redshift_kwargs):
    #method='Batten2020', cosmology="Planck18", sigma=1,
    #                           scatter_percentage=0, subtract_host=False, lookup_table=None):
        """
        Calculates the mean redshift and the confidence interval of
        an FRB from its :attr:`dm`, :attr:`dm_excess` or
        :attr:`dm_excess` - :attr:`dm_host_est`.

        Parameters
        ----------
        method : str, optional
            The dispersion meausre-redshift relation to use when
            calculating the redshift. Avaliable methods:  %(meth)s.
            Default: 'Batten2020'

        cosmology : str, optional
            The cosmology to assume when calculating the redshift.
            This value is overided if using a hydrodynamic method.
            Avaliable cosmologies: %(cosmo)s. Default: 'Planck18'

        sigma : int (1, 2, 3, 4, 5), optional
            The width of the confidence interval in units of standard
            deviation. `sigma=1` is the 68\% confidence interval.

        scatter_percentage : float, optional
            The amount of line of
            Default: 0

        subtract_host : bool, optional
            Subtract :attr:`dm_host_est` from the :attr:`dm_excess`
            before calculating the redshift. This is is used to account
            for the dispersion measure that arises from the FRB host
            galaxy. Default: False

        lookup_table : str or None, optional
            The path to the lookup table file. If ``lookup_table=None``
            a table will attempted to be loaded from the data directory
            based on the method name. Default: *None*

        Returns
        -------
        float
            The extimated redshift of the FRB.

        float
            The redshift confidence interval the FRB.

        Notes
        -----
        The methods_ section in the documentation has a description for
        each of the builtin methods.

        The cosmology_ section in the documentation has a list of the
        cosmological parameters for each cosmology

        .. _cosmology: https://fruitbat.readthedocs.io/en/latest/user_guide/method_and_cosmology.html#cosmology
        .. _methods: https://fruitbat.readthedocs.io/en/latest/user_guide/method_and_cosmology.html#methods

        """
        redshift = self.calc_redshift(**calc_redshift_kwargs)


        if method in methods.methods_hydrodynamic():
            if method == "Batten2021":
                zvals, pdf, dz = self.calc_redshift_pdf(method="Batten2021")

                conf_lower_lim, conf_upper_lim = utils.sigma_to_pdf_percentiles(sigma=sigma)

                conf_int_lower = utils.calc_z_from_pdf_percentile(zvals, pdf, percentile=conf_lower_lim)
                conf_int_upper = utils.calc_z_from_pdf_percentile(zvals, pdf, percentile=conf_upper_lim)

        self.z_conf_int_lower = conf_int_lower
        self.z_conf_int_upper = conf_int_upper

        return self.z, (self.z_conf_int_lower, self.z_conf_int_upper)


    def calc_skycoords(self, frame=None):
        """
        Calculates the skycoord position on the sky of the FRB from
        (:attr:`raj`, :attr:`decj`) or (:attr:`gl`, :attr:`gb`). If
        both (:attr:`raj`, :attr:`decj`) and (:attr:`gl`, :attr:`gb`)
        are given, preference is given to (:attr:`raj`, :attr:`decj`).

        Parameters
        ----------
        frame: str or None, optional
            The type of coordinate frame. If ``frame = None`` then
            :meth:`calc_skycoords` will use a default frame based on
            the coordinates given. If :attr:`raj` and :attr:`decj` are
            given the default frame is 'icrs'. If :attr:`gl` and
            :attr:`gb` are given the default frame is 'galactic'.
            Default: *None*

        Returns
        -------
        :obj:`astropy.coordinates.SkyCoord`
            The sky coordinates of the FRB.
        """

        if self.raj is not None and self.decj is not None:
            if frame is None:
                frame = "icrs"
            skycoords = SkyCoord(ra=self.raj, dec=self.decj, frame=frame,
                                 unit=(u.hourangle, u.deg))

        elif self.gl is not None and self.gb is not None:
            if frame is None:
                frame = "galactic"
            skycoords = SkyCoord(self.gl, self.gb, frame=frame,
                                 unit=u.deg)
        else:
            raise ValueError("To calculate skycoords either (raj and decj)"
                             "or (gl, gb) must be provided")

        return skycoords

    def calc_dm_excess(self):
        """
        Calculates the dispersion measure excess of the FRB by
        subtracting the DM contribution from the Milky Way.

        Returns
        -------
        :obj:`astropy.units.Quantity`
            The dispersion measure excess.

        Notes
        -----
        :math:`\\rm{DM_{excess}}` is calculated as follows:

        .. math::

            \\rm{DM_{excess} = DM - DM_{galaxy}}
        """
        dm_excess = self.dm.value - self.dm_galaxy.value
        if dm_excess < 0:
            print("dm_excess < 0: This implies that the DM estimate "
                  "from the Milky Way is higher than the observed DM. "
                  "Setting dm_excess = 0")
            self.dm_excess = 0
        else:
            self.dm_excess = dm_excess
        return dm_excess

    def calc_dm_galaxy(self, model='ymw16', include_halo=False, return_tau_sc=False):
        """
        Calculates the dispersion measure contribution of the Milky Way
        from either (:attr:`raj`, :attr:`decj`) or (:attr:`gl`,
        :attr:`gb`). Uses the YMW16 model of the Milky Way free
        electron column density.

        Parameters
        ----------
        model : 'ymw16' or 'ne2001', optional
            The Milky Way dispersion measure model. To use 'ne2001' you
            will need to install the python port. See
            https://fruitbat.readthedocs.io/en/latest/user_guide/ne2001_installation.html
            Default: 'ymw16'

        include_halo : bool, optional
            Include the DM of the galactic halo which isn't included
            in NE2001 or YMW16. This used the YT2020 halo model. Default: False

        return_tau_sc : bool, optional
            Return the scattering timescale in addition to the DM.
            Default: False

        Returns
        -------
        dm_galaxy: :obj:`astropy.units.Quantity`
            The dispersion measure contribution from the Milky Way of
            the FRB.

        tau_sc: :obj:`astropy.units.Quantity`, optional
            The scattering timescale at 1 GHz (s). Only returns if
            :attr:`return_tau_sc` is `True`.

        """
        YMW16_options = ["ymw16", "YMW16"]
        NE2001_options = ["ne2001", "NE2001"]

        if model in YMW16_options:
            model = 'YMW16'
        elif model in NE2001_options:
            model = 'NE2001'
        else:
            raise ValueError("'{}' is not a valid galactic DM model".format(model))

        # Check to make sure some of the keyword are not None
        coord_list = [self.skycoords, self.raj, self.decj, self.gl, self.gb]
        if all(val is None for val in coord_list):
            raise ValueError("""Can not calculate dm_galaxy since coordinates
                             for FRB burst were not provided. Please provide
                             (raj, decj) or (gl, gb) coordinates.""")

        # Calculate skycoords position if it
        elif (self.skycoords is None and
                (self.raj is not None and self.decj is not None) or
                (self.gl is not None and self.gb is not None)):

            self._skycoords = self.calc_skycoords()

        dm_galaxy, tau_sc = pygedm.dist_to_dm(
            gl=self._skycoords.galactic.l,
            gb=self._skycoords.galactic.b,
            dist=25000.0,
            method=model)

        if include_halo:
            self.galaxy_halo = pygedm.calculate_halo_dm(
                gl=self._skycoords.galactic.l,
                gb=self._skycoords.galactic.b,
                method='yt2020'
            )

        else:
            self.galaxy_halo = 0 * u.pc * u.cm**(-3)

        self.dm_galaxy_model = model
        self.dm_galaxy = dm_galaxy.value + self.galaxy_halo.value
        self.tau_sc = tau_sc.value
        self.calc_dm_excess()

        if return_tau_sc:
            return self.dm_galaxy, self.tau_sc
        else:
            return self.dm_galaxy


    def calc_dm_igm(self):
        """
        Calculates the dispersion measure of the intergalactic medium
        along the line-of-sight of the FRB. This can only be done if
        the redshift and dispersion measure contribution of the FRB
        host galaxy is known.

        Returns
        -------
        :obj:`astropy.units.Quantity`
            The dispersion measure contribution of the IGM.

        Notes
        -----
        :math:`\\rm{DM_{IGM}}` is calculated as follows:

        .. math::

           \\rm{DM_{IGM}} = \\rm{DM_{excess}} -
           \\frac{\\rm{DM_{host, loc}}}{1 + z}
        """

        if self.z_host is None:
            err_msg = ("z_host is None. Provide a non zero value for the "
                       "FRB host redshift")
            raise ValueError(err_msg)

        if self.dm_host_loc == 0.0:
            err_msg = ("dm_host_loc = 0. The dm_igm will be the same as "
                       "dm_excess. Provide a non-zero value for dm_host_loc")
            raise ValueError(err_msg)

        dm_igm = self.dm_excess - (self.dm_host_loc / (1 + self.z_host))
        self.dm_igm = dm_igm.value
        return self.dm_igm

    def calc_fluence(self):
        """
        Calculates the observed fluence of the FRB. This requires
        :attr:`width` and :attr:`peak_flux` to not be *None*.

        Returns
        -------
        :obj:`astropy.units.Quantity` or None:
            The fluence of the FRB.

        Notes
        -----
        The fluence (:math:`\\rm{F_{obs}}`) is calculated as follows:

        .. math::

            F_{obs} = W_{obs} S_{\\nu, p}

        Where :math:`W_{obs}` is the with of the FRB pulse and
        :math:`S_{\\nu, p}` is the specific peak flux.
        """
        if (self.width is None) or (self.peak_flux is None):
            err_msg = ("calc_fluence requires both width and peak_flux "
                       "to not be None")
            raise ValueError(err_msg)

        fluence = self.width * self.peak_flux
        self.fluence = fluence.value
        return fluence

    def calc_luminosity_distance(self):
        """
        Calculates the luminosity distance to the FRB based on its
        redshift. To calculate a luminosity distance either call
        :meth:`calc_redshift` first to determine a redshift estimate
        or provide a value for :attr:`dm_host`.

        Returns
        -------
        :obj:`astropy.units.Quantity`
            The luminosity distance to the FRB.
        """

        if self.z is None:
            raise ValueError(
                """ Can not calculate luminosity distance without a redshift.
                Use calc_redshift() to calculate the FRB redshift or provide
                a value for z_host.
                """)

        cosmo = cosmologies.cosmology_functions()[self.cosmology]
        return cosmo.luminosity_distance(self.z)

    def calc_comoving_distance(self):
        """
        Calculates the comoving distance to the FRB based on its
        redshift. To calculate a comoving distance either call
        :meth:`calc_redshift` first to determine a redshift estimate
        or provide a value for :attr:`dm_host`.

        Returns
        -------
        :obj:`astropy.units.Quantity`
            The comoving distance to the FRB.
        """
        if self.z is not None:
            z_sample = self.z
        elif self.z_host is not None:
            z_sample = self.z_host
        else:
            raise ValueError(
                """ Can not calculate comoving distance without a redshift.
                Use calc_redshift() to calculate the FRB redshift or provide
                a value for z_host.
                """)

        cosmo = cosmologies.cosmology_functions()[self.cosmology]
        return cosmo.comoving_distance(z_sample)

    def calc_luminosity(self, use_bandwidth=False):
        """
        Calculates the isotropic peak luminosity of the FRB. This is
        the upper limit to the the true peak luminosity since the
        luminosity distance required is also an upper limit to the true
        luminosity distance.

        Parameters
        ----------
        use_bandwidth: bool, optional
            The default method of calculating the luminosity of a FRB
            uses :attr:`obs_freq_central` as described in Zhang 2018.
            However some estimates of the FRB luminosity instead use
            :attr:`obs_bandwidth` (see Law et al. 2017). Set to
            ``True`` to use :attr:`obs_bandwidth` instead of
            :attr:`obs_freq_central`. Default: False

        Returns
        -------
        :obj:`astropy.units.Quantity`
            The estimated isotropic peak luminosity of the FRB in units
            of ergs/s.

        Notes
        -----
        The luminosity of a FRB is calculated following that in
        Zhang 2018.

        .. math::

            L_{FRB} \\simeq 4 \\pi D_{\\rm{L}}^2 S_{\\nu, p} \\nu_c

        Where :math:`S_{\\nu, p}` is the specific peak flux of the FRB,
        :math:`\\nu_c` is the central observing frequency,
        :math:`D_{\\rm{L}}` is the luminosity distance.

        In Law et al. 2017, the bandwidth (:math:`B`) is used instead
        of :math:`\\nu_c` to estimate luminosity.

        """
        D = self.calc_luminosity_distance()

        if self.peak_flux is None:
            err_msg = ("Can not calculate energy without peak_flux. Provide "
                       "peak_flux before calculating luminosity.")
            raise ValueError(err_msg)
        else:
            S = self.peak_flux

        if use_bandwidth:
            if self.obs_bandwidth is None:
                err_msg = ("Can not calculate energy without observing "
                           "bandwidth. Provide obs_bandwidth before "
                           "calculating energy.")
                raise ValueError(err_msg)
            B = self.obs_bandwidth
            lum = 4 * np.pi * D**2 * S * B
        else:
            if self.obs_freq_central is None:
                err_msg = ("Can not calculate energy without observing "
                           "frequency. Provide obs_freq_central before "
                           "calculating energy.")
                raise ValueError(err_msg)
            nu = self.obs_freq_central
            lum = 4 * np.pi * D**2 * S * nu

        return lum.to('erg s**-1')

    def calc_energy(self, use_bandwidth=False):
        """
        Calculates the isotropic energy of the FRB. This is the upper
        limit to the the true energy since the luminosity distance
        required is also an upper limit to the true luminosity
        distance.

        Parameters
        ----------
        use_bandwidth: bool, optional
            The default method of calculating the energy of a FRB uses
            :attr:`obs_freq_central` as described in Zhang 2018.
            However some estimates of the FRB energy instead use
            :attr:`obs_bandwidth` (see Law et al. 2017). Set to
            ``True`` to use :attr:`obs_bandwidth` instead of
            :attr:`obs_freq_central`. Default: False

        Returns
        -------
        :obj:`astropy.units.Quantity`
            The estimated isotropic energy of the FRB, in units of
            ergs.

        Notes
        -----
        The energy of a FRB is calculated following that in Zhang 2018.

        .. math::

            E_{FRB} \\simeq \\frac{4 \\pi}{1 + z} D_{\\rm{L}}^2 F_{obs} \\nu_c

        Where :math:`F_{obs}` is the fluence of the FRB, :math:`\\nu_c`
        is the central observing frequency, :math:`D_{\\rm{L}}` is the
        luminosity distance and :math:`z` is the estimated redshift.

        In Law et al. 2017, the bandwidth (:math:`B`) is used instead
        of :math:`\\nu_c` to estimate energy.
        """

        D = self.calc_luminosity_distance()

        if self.fluence is None:
            err_msg = ("Can not calculate energy without fluence. Provide "
                       "fluence or peak_flux and width before calculating "
                       "energy.")
            raise ValueError(err_msg)
        else:
            fluence = self.fluence

        if use_bandwidth:
            if self.obs_bandwidth is None:
                err_msg = ("Can not calculate energy without observing "
                           "bandwidth. Provide obs_bandwidth before "
                           "calculating energy.")
                raise ValueError(err_msg)
            else:
                bandwidth = self.obs_bandwidth
            energy = fluence * bandwidth * 4 * np.pi * D**2 * (1 + self.z)**-1

        else:
            if self.obs_freq_central is None:
                err_msg = ("Can not calculate energy without observing "
                           "frequency. Provide obs_freq_central before "
                           "calculating energy.")
                raise ValueError(err_msg)
            else:
                nu = self.obs_freq_central
            energy = fluence * nu * 4 * np.pi * D**2 * (1 + self.z)**-1

        return energy.to("erg")

    def _set_value_units(self, value, unit=None, non_negative=False):
        """
        Converts ``value`` into a :obj:`~astropy.unit.Quantity` with
        units. Also checks if ``value`` has existing units and the
        units are convertable with ``unit``. If ``value`` is a
        :obj:`~astropy.unit.Quantity` with convertable units then
        ``value`` is returned unchanged.

        Parameters
        ----------
        value: float, int, or :obj:`~astropy.unit.Quantity`
            The input value.

        unit: :obj:`astropy.unit` or None
            The unit to set for ``value``. If ``unit = None`` then
            ``value`` will become a dimensionless quantity.

        non_negative: bool, optional
            Raise an error if the value is negative. This is used to
            verify values that can only be positive are correctly
            specified. Default: *False*

        Return
        ------
        var : :obj:`~astropy.unit.Quantity` or None
            The output quantity with units. If ``value=None``
            returns *None*
        """
        # Check if the value already has units
        if not isinstance(value, u.Quantity):
            if value is None:
                var = None

            elif non_negative and value < 0.0:
                raise ValueError("Value must be greater than zero.")

            elif unit is None:
                var = value * u.dimensionless_unscaled

            else:
                var = value * unit
        else:
            # If passed units check that they are the correct units by checking
            # that they are equivalent

            if value.unit.is_equivalent(unit):
                var = value

            else:
                err_msg = ("Quantity expected to have units of {} instead was "
                           "passed units {} which is not convertable".format(
                            unit, value.unit))
                raise ValueError(err_msg)

        return var

    @property
    def name(self):
        """
        name: str
            The name of the FRB object.
        """
        return self._name

    @name.setter
    def name(self, value):
        if isinstance(value, str):
            self._name = value
        elif value is None:
            self._name = None
        else:
            self._name = str(value)

    @property
    def dm(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The observed dispersion measure of the FRB. This is without
            Milky Way or host galaxy subtraction.
        """
        return self._dm

    @dm.setter
    def dm(self, value):

        utils.check_type("dm", value, str, desire=False)


        self._dm = self._set_value_units(value, unit=u.pc * u.cm**-3,
                                         non_negative=True)

    @property
    def dm_galaxy(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The galactic component of the dispersion measure. This is
            quantity is calculated using :meth:`calc_dm_galaxy`.
        """
        return self._dm_galaxy

    @dm_galaxy.setter
    def dm_galaxy(self, value):
        self._dm_galaxy = self._set_value_units(value, unit=u.pc * u.cm**-3,
                                                non_negative=True)

    @property
    def dm_excess(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The dispersion measure with the Milky Way component
            subtracted. This value approximates the dispersion measure
            of the intergalatic medium by assuming that the host galaxy
            contributes zero to the observed dispersion measure. This
            quantity is calculated using :meth:`calc_dm_excess`.
        """
        return self._dm_excess

    @dm_excess.setter
    def dm_excess(self, value):
        self._dm_excess = self._set_value_units(value, u.pc * u.cm**-3,
                                                non_negative=True)


    @property
    def dm_host_est(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The estimated dispersion measure from the FRB host galaxy.
        """
        return self._dm_host_est

    @dm_host_est.setter
    def dm_host_est(self, value):
        self._dm_host_est = self._set_value_units(value, u.pc * u.cm**-3,
                                                  non_negative=True)

    @property
    def dm_host_loc(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The dispersion measure from a localised FRB host galaxy.
        """
        return u.Quantity(self._dm_host_loc, u.pc * u.cm**-3)

    @dm_host_loc.setter
    def dm_host_loc(self, value):
        self._dm_host_loc = self._set_value_units(value, u.pc * u.cm**-3,
                                                  non_negative=True)

    @property
    def dm_index(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The dispersion measure index of the burst.
        """
        return self._dm_index

    @dm_index.setter
    def dm_index(self, value):
        self._dm_index = self._set_value_units(value)

    @property
    def z(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The estimated redshift of the burst. By default this
            assumes that the entire :attr:`dm_excess` arrives from the
            IGM and the host galaxy of the FRB and any surrounding
            material contribute nothing to the total DM. This should be
            taken as an upper limit to the bursts true redshift. To
            provide an estimate of the DM contribution due to he host
            galaxy, set :attr:`dm_host_est` to a non-zero value and use
            ``subract_host=True`` when calling :meth:`calc_redshift()`.
        """
        return self._z

    @z.setter
    def z(self, value):
        self._z = self._set_value_units(value, non_negative=True)


    @property
    def z_conf_int_lower(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The lower bound of the redshift confidence interval.
        """
        return self._z_conf_int_lower

    @z_conf_int_lower.setter
    def z_conf_int_lower(self, value):
        self._z_conf_int_lower = self._set_value_units(value, non_negative=True)


    @property
    def z_conf_int_upper(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The upper bound of the redshift confidence interval.
        """
        return self._z_conf_int_upper

    @z_conf_int_upper.setter
    def z_conf_int_upper(self, value):
        self._z_conf_int_upper = self._set_value_units(value, non_negative=True)


    @property
    def z_host(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The redshift of the localised FRB host galaxy. Note that
            this an observed quantity, not the estimated redshift
            :attr:`z` calculated with :meth:`calc_redshift()`.
        """
        return self._z_host

    @z_host.setter
    def z_host(self, value):
        self._z_host = self._set_value_units(value)

    @property
    def scatt_index(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The scattering index of the burst.
        """
        return self._scatt_index

    @scatt_index.setter
    def scatt_index(self, value):
        self._scatt_index = self._set_value_units(value)

    @property
    def width(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The observed width of the pulse.
        """
        return self._width

    @width.setter
    def width(self, value):
        self._width = self._set_value_units(value, u.ms, non_negative=True)

    @property
    def peak_flux(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The observed peak flux density of the FRB.
        """
        return self._peak_flux

    @peak_flux.setter
    def peak_flux(self, value):
        self._peak_flux = self._set_value_units(value, u.Jy, non_negative=True)

    @property
    def fluence(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The observed fluence of the FRB. This is calculated from
            :meth:`calc_fluence`.
        """
        return self._fluence

    @fluence.setter
    def fluence(self, value):
        self._fluence = self._set_value_units(value, u.Jy * u.ms,
                                              non_negative=True)

    @property
    def obs_bandwidth(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The observing bandwidth of the FRB.
        """
        return self._obs_bandwidth

    @obs_bandwidth.setter
    def obs_bandwidth(self, value):
        self._obs_bandwidth = self._set_value_units(value, u.MHz,
                                                    non_negative=True)

    @property
    def obs_freq_central(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The central observing frequency of the FRB.
        """
        return self._obs_freq_central

    @obs_freq_central.setter
    def obs_freq_central(self, value):
        self._obs_freq_central = self._set_value_units(value, u.MHz,
                                                       non_negative=True)

    @property
    def raj(self):
        """
        :obj:`astropy.coordinates.Longitude` or None:
            The right accension in J2000 coordinates of the best
            estimate of the FRB position. This is given in the IRCS
            frame.
        """
        return self._raj

    @raj.setter
    def raj(self, value):
        if value is None:
            self._raj = None
        else:
            self._raj = coord.angles.Longitude(value, u.hourangle)

    @property
    def decj(self):
        """
        :obj:`astropy.coordinates.Latitude` or None:
            The declination in J2000 coordinates of the best estimate
            of the FRB position. This is given in the IRCS frame.
        """
        return self._decj

    @decj.setter
    def decj(self, value):
        if value is None:
            self._decj = None
        else:
            self._decj = coord.angles.Latitude(value, u.deg)

    @property
    def gl(self):
        """
        :obj:`astropy.coordinates.Longitude` or None:
            The longitude in galactic coordinates of the best estimate
            of the FRB position.
        """
        return self._gl

    @gl.setter
    def gl(self, value):
        if value is None:
            self._gl = None
        else:
            self._gl = coord.angles.Longitude(value, u.deg)

    @property
    def gb(self):
        """
        :obj:`astropy.coordinates.Latitude` or None:
            The latitude in galactic coordinates of the best estimate
            of the FRB position.
        """
        return self._gb

    @gb.setter
    def gb(self, value):
        if value is None:
            self._gb = None
        else:
            self._gb = coord.angles.Latitude(value, u.deg)

    @property
    def skycoords(self):
        """
        :obj:`astropy.coordinates.SkyCoord` or None:
            The skycoords of the FRB. This is calculated from either
            (:attr:`raj`, :attr:`decj`) or (:attr:`gl`, :attr:`gb`)
            using :meth:`calc_skycoords`.
        """
        return self._skycoords

    @skycoords.setter
    def skycoords(self, value):
        if value is None:
            self._skycoords = None
        else:
            self._skycoords = SkyCoord(value)

    @property
    def dm_igm(self):
        """
        :obj:`astropy.units.Quantity` or None:
            The estimated disperison measure from the IGM. This can be
            calculated for a FRB for which the redshift of its host
            galaxy is known. This value is determined using
            :meth:`calc_dm_igm`.
        """
        return self._dm_igm

    @dm_igm.setter
    def dm_igm(self, value):
        self._dm_igm = self._set_value_units(value, u.pc * u.cm**-3,
                                             non_negative=True)

    @property
    def utc(self):
        """
        :obj:`astropy.time.Time` or None
        The UTC time of the burst.
        """
        return self._utc

    @utc.setter
    def utc(self, value):
        if value is None:
            self._utc = None
        else:
            self._utc = Time(value, format='isot', scale='utc')

    @property
    def snr(self):
        """
        :obj:`astropy.units.Quantity`:
        The signal-to-noise  ratio of the burst.
        """
        return self._snr

    @snr.setter
    def snr(self, value):
        self._snr = self._set_value_units(value, non_negative=True)

    @property
    def cosmology(self):
        """
        str or None:
            The cosmology used to calculate redshift.
        """
        return self._cosmology

    @cosmology.setter
    def cosmology(self, value):
        if isinstance(value, str):
            self._cosmology = value
        else:
            self._cosmology = str(value)

    @property
    def method(self):
        """
        str or None:
            The method used to calculate redshift.
        """
        return self._method

    @method.setter
    def method(self, value):
        if isinstance(value, str):
            self._method = value
        else:
            self._method = str(value)


    @property
    def dm_galaxy_model(self):
        """
        str or None:
            The method used to calculate redshift.
        """
        return self._dm_galaxy_model

    @dm_galaxy_model.setter
    def dm_galaxy_model(self, value):
        if isinstance(value, str):
            self._dm_galaxy_model = value
        else:
            self._dm_galaxy_model = str(value)