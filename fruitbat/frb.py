"""
Provides the definition of the Frb class and methods.
"""
from e13tools import docstring_substitute

from astropy.coordinates import SkyCoord
import astropy.units as u

import pyymw16 as ymw16

from . import estimate
from . import cosmology
from ._fruitbatstrings import (dm_units_doc)

__all__ = ["Frb"]

@docstring_substitute(dm_units=dm_units_doc)
class Frb(object):
    """
    Defines the :class:`~Frb` class in the **fruitbat** package.

    Arguments
    ---------
    name : str
        The name of the Frb object.

    dm : float
        The observed dispersion measure of the FRB. This is without Milky Way
        or host galaxy subtraction.
        Units: %(dm_units)s

    Keyword Arguments
    -----------------
    raj : str or None, optional
        The right ascension in J2000 coordinates of the best estimate of the
        FRB position. Default: *None*
        
    decj : str or None, optional
        The declination in J2000 coordinates of the best estimate of the
        FRB position. Default: *None*        

    gl : str or None, optional
        The Galactic longitude in degrees of the best estimate of the FRB 
        position. Default: *None*
    
    gb : str or None, optional
        The Galactic latitude in degrees of the best estimate of the FRB 
        position. Default: *None*
  
    Other Parameters
    ----------------
    dm_galaxy : float, optional
        The modelled contribution to the FRB DM by electrons in the Milky Way.
        Units: %(dm_units)s Default: 0.0

    dm_excess : float or None, optional
        The DM excess of the FRB over the estimated Galactic DM. If
        :attr:`dm_excess` is *None*, then :attr:`dm_excess` is calculated 
        automatically with :meth:`calc_dm_excess()`.
        Units: %(dm_units)s Default: *None*

    z_host : float or None, optional
        The observed redshift of the localised FRB host galaxy.
        Default: *None*

    dm_host_est : float, optional
        The estimated contribution to the measured FRB DM from originating from
        the FRB's host galaxy. This value is the amount of DM the host galaxy
        contributes to the observed DM, *not* the DM of the host galaxy.
        Units: %(dm_units)s Default: 0.0

    dm_host_loc : float, optional
        The dispersion measure of a localised FRB host galaxy. This value is
        *not* the contribution to the observed DM, but the DM at the host 
        galaxy. The observed DM is :attr:`dm_host_loc` but attenuated by a 
        factor of (1 + z).
        Units: %(dm_units)s Default: 0.0

    dm_index : float or None, optional
        The dispersion measure index of the burst :math:`\\alpha` such that
        :math:`\\rm{DM} \\propto \\nu^{-\\alpha}` Default: *None*

    scatt_index : float or None, optional
        The scattering index (:math:`\\beta`) of the FRB pulse. The
        scattering index describes how the width (:math:`\\rm{W}`) of the
        FRB pulse evolves with frequency :math:`\\nu` such that
        :math:`\\rm{W} \\propto \\nu^{-\\beta}`. Default: *None*

    snr : float or None, optional
        The signal-to-noise of the burst.
        Default: *None*

    w_obs : float or None, optional
        The observed width of the pulse obtained by a pulse fitting algorithm.
        Units: :math:`\\rm{ms}` Default: *None*

    s_peak_obs : float or None, optional
        The observed peak flux density of the burst.
        Units: :math:`\\rm{Jy}` Default: *None*

    f_obs : float or None, optional
        The observed fluence of the FRB. If :attr:`f_obs` is *None* and both
        :attr:`w_obs` and :attr:`s_peak_obs` are not *None* then :attr:`f_obs` 
        is automatically calculated by :attr:`w_obs` x :attr:`s_peak_obs`
        Units: :math:`\\rm{Jy\\ ms}` Default: *None*

    utc : str or None, optional
        The UTC time of the FRB Burst. Default: *None*

    dm_uncert : float, optional
        The uncertainty in the dispersion measure.
        Units: %(dm_units)s Default: 0.0

    z_uncert : float, optional
        The uncertainty in the redshift of the FRB. Default: 0.0

    """

    def __init__(self, name, dm, *, raj=None, decj=None, gl=None, gb=None,
                 dm_galaxy=0.0, dm_excess=None, z_host=None, dm_host_est=0.0,
                 dm_host_loc=0.0, dm_index=None, scatt_index=None, snr=None, 
                 w_obs=None, s_peak_obs=None, f_obs=None, utc=None,
                 dm_uncert=0.0, z_uncert=0.0):

        # TO DO:
        # There are a few other parameters that I should add:
        # dm_index


        self._name = name

        self._dm = float(dm)
        if self._dm < 0.0:
            raise ValueError("Dispersion Measure can not be negative.")

        self._dm_uncert = float(dm_uncert)
        if self._dm_uncert < 0.0:
            raise ValueError("Dispersion Measure can not be negative.")

        self._dm_galaxy = float(dm_galaxy)
        if self._dm_galaxy < 0.0:
            raise ValueError("Dispersion Measure can not be negative.")

        self._dm_host_est = float(dm_host_est)
        if self._dm_host_est < 0.0:
            raise ValueError("Dispersion Measure can not be negative.")

        self._dm_host_loc = float(dm_host_loc)
        if self._dm_host_loc < 0.0:
            raise ValueError("Dispersion Measure can not be negative.")

        # Calculate dm_excess from existing parameters if it is not given.
        if not dm_excess:
            self.calc_dm_excess()
        else:
            self._dm_excess = dm_excess

        self._dm_index = dm_index
        self._z_host = z_host
        self._z_uncert = z_uncert
        self._scatt_index = scatt_index
        self._snr = snr
        self._w_obs = w_obs
        self._s_peak_obs = s_peak_obs
        self._utc = utc
        self._raj = raj
        self._decj = decj
        self._gl = gl
        self._gb = gb
        self._z = None

        # Calculate F_obs if s_peak_obs and w_obs are given
        if (not f_obs) and (s_peak_obs and w_obs):
            self.calc_f_obs()
        else:
            self._f_obs = f_obs

        # Calculate the skycoords of the FRB from (raj, decj) or (gl, gb)
        if (raj and decj) is not None or (gl and gb) is not None:
            self._skycoords = self.calc_skycoords()
            self._raj = self.skycoords.icrs.ra
            self._decj = self.skycoords.icrs.dec
            self._gl = self.skycoords.galactic.l
            self._gb = self.skycoords.galactic.b
        else:
            self._skycoords = None

    def __repr__(self):
        return 'Frb({0})'.format(vars(self))

    @docstring_substitute(methods=estimate.methods(string=True), 
                          cosmo=cosmology.keys(string=True))
    def calc_redshift(self, method='inoue2004', cosmology="planck2018",
        subtract_host=False):
        """
        Calculate the redshift of the FRB from its :attr:`dm`, :attr:`dm_excess` 
        or :attr:`dm_excess` - :attr:`dm_host_est`.

        Parameters
        ----------
        method : str, optional
            The approximation to use when calculating the redshift.
            Avaliable methods:  %(methods)s

        cosmology : str, optional
            The method `inoue2004` has the option to choose which cosmology
            to assume when performing the redshift estimation.
            Avaliable cosmologies: %(cosmo)s

        subtract_host : bool, optional
            Subtract :attr:`dm_host_est` from the :attr:`dm_excess` before 
            calculating the redshift. This is is used to account for the 
            dispersion measure that arises from the FRB host galaxy.

        Returns
        -------
        float
            The redshift of the FRB.

        Notes
        -----
        The methods_ section in the documentation has a discription of each 
        methods and where they should apply.

        The cosmology_ section of the documentation has a list of the 
        cosmological parameters used in each cosmology method.

        .. _cosmology: https://fruitbat.readthedocs.io/en/latest/user_guide/method_and_cosmology.html#cosmology
        .. _methods: https://fruitbat.readthedocs.io/en/latest/user_guide/method_and_cosmology.html#methods
        """
        if not isinstance(subtract_host, bool):
            raise ValueError("subtract_host must be of type bool.")

        if subtract_host:
            input_dm = self.dm_excess - self.dm_host_est
        else:
            input_dm = self.dm_excess

        z  = estimate.redshift(dm=input_dm, method=method, 
                               cosmology=cosmology)
        self._z = z
        return z


    def calc_skycoords(self):
        """
        Calculates the skycoord position on the sky of the FRB from 
        (:attr:`raj`, :attr:`decj`) or (:attr:`gl`, :attr:`gb`).

        Returns
        -------
        astropy.coordinates.sky_coordinate.SkyCoord
            The sky coordinates of the FRB.
        """

        if self.raj is not None and self.decj is not None:
            skycoords = SkyCoord(self.raj, self.decj, frame="icrs", 
                                 unit=(u.hourangle, u.deg))

        elif self.gl is not None and self.gb is not None:
            skycoords = SkyCoord(self.gl, self.gb, frame="galactic", 
                                unit=u.deg)
        else:
            raise ValueError("To calculate skycoords either (raj and decj)"
                             "or (gl, gb) must be provided")


        return skycoords


    def calc_dm_excess(self):
        """
        Calculates the dispersion measure excess of the FRB by subtracting
        the DM contribution from the Milky Way.

        Returns
        -------
        dm_excess : float
            The dispersion measure excess.

        Notes
        -----
        :math:`\\rm{DM_{excess}}` is calculated as follows:

        .. math::

            DM_{excess} = DM - DM_{galaxy}
        """
        dm_excess = self.dm - self.dm_galaxy
        self._dm_excess = dm_excess
        return dm_excess


    def calc_dm_galaxy(self, model='ymw16'):
        """
        Calculates the dispersion measure contribution of the Milky Way from
        either (:attr:`raj`, :attr:`decj`) or (:attr:`gl`, :attr:`gb`).

        Parameters
        ----------
        model : str, optional
            The Milky Way dispersion measure model. Default: ymw16
        
        """

        # Since the YMW16 model only gives you a dispersion measure out to a 
        # distance within the galaxy, to get the entire DM contribution of the 
        # galaxy we need to specify the furthest distance in the YMW16 model.
        max_galaxy_dist = 25000  # units: pc

        if not self.skycoords and (self.raj and self.decj) or (self.gl, self.gb):
            self._skycoords = self.calc_skycoords()
        elif not self.skycoords and not (self.raj and self.decj) and not (self.gl and self.gb):
            raise ValueError("""Can not calculate dm_galaxy since coordinates
                             for FRB burst were not provided. Please provide
                             (raj, decj) or (gl, gb) coordinates.""")

        dm_galaxy, tau_sc = ymw16.dist_to_dm(self.skycoords.galactic.l, 
                                             self.skycoords.galactic.b, 
                                             max_galaxy_dist)

        self._tau_sc = tau_sc
        self._dm_galaxy = dm_galaxy.value
        dm_excess = self.calc_dm_excess()
        return dm_galaxy.value


    def calc_dm_igm(self):
        """
        Calculates the dispersion measure of the intergalactic medium along the
        line-of-sight of the FRB. This can only be done if the redshift and 
        dispersion measure contribution of the FRB host galaxy is known.
        
        Returns
        -------
        dm_igm : float
            The dispersion measure contribution of the IGM.

        Notes
        -----
        :math:`DM_{IGM}` is calculated as follows:

        .. math::

           DM_{IGM} = DM_{excess} - \\frac{DM_{host,loc}}{1 + z} 
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
        self._dm_igm = dm_igm
        return dm_igm


    def calc_f_obs(self):
        """
        Calculates the observed fluence of the FRB. This requires :attr:`w_obs`
        and :attr:`s_peak_obs` to not be *None*.

        Returns
        -------
        float
            The fluence of the FRB.

        Notes
        -----
        :math:`\\rm{F_{obs}}` is calculated as follows:

        .. math::

            \\rm{F_{obs} = W_{obs} \\times S_{peak, obs}}
        """

        if (not self.w_obs) or (not self.s_peak_obs):
            err_msg = ("calc_f_obs requires both w_obs and s_peak_obs "
                       "to not be None")
            raise ValueError(err_msg)

        f_obs = self.w_obs * self.s_peak_obs
        self._f_obs = f_obs
        return f_obs


#        def calc_d_comov(self, cosmology):



#    def pulse_width(self, freq):
#        """
#        The width of the pulse at a given frquency using the scattering index.
#
#        Parameters
#        ----------
#        freq: float
#            The frequency at
#
#        Returns
#        -------
#        float:
#            The width of the pulse. 
#        """
#        coeff = 
#        return coeff * freq**(-self.scatt_index)


    @property
    def name(self):
        """
        str: 
            The name of the FRB object.
        """
        return self._name   

    @property
    def dm(self):
        """
        float: 
            The observed dispersion measure of the FRB.
        """
        return self._dm    
 
#    @property
#    def dm_uncert(self):
#        """
#        float: The uncertainty in the observed dispersion measure of the FRB.
#        """
#        return self._dm_uncert

    @property
    def dm_galaxy(self):
        """
        float: 
            The Milky Way component of the dispersion measure.
        """
        return self._dm_galaxy    
    
    @property
    def dm_excess(self):
        """
        float: 
            The dispersion measure with the Milky Way component subtracted.
        """
        return self._dm_excess
    
    @property
    def dm_host_est(self):
        """
        float: 
            The dispersion measure from the FRB host galaxy
        """
        return self._dm_host_est

    @property
    def dm_host_loc(self):
        """
        float: 
            The dispersion measure from a localised FRB host galaxy
        """
        return self._dm_host_loc

#    @property
#    def dm_index(self):
#        """
#        float: The dispersion measure index of the burst.
#        """
#        return self._dm_index

    @property
    def z(self):
        """    
        float or None:
            The estimated redshift of the burst. By default this assumes that 
            the entire :attr:`dm_excess` arrives from the IGM and the host 
            galaxy of the FRB and any surrounding material contribute nothing 
            to the total DM. This should be taken as an upper limit to the 
            bursts true redshift. To provide an estimate of the DM contribution
            due to he host galaxy, set :attr:`dm_host_est` to a non-zero value 
            and use ``subract_host=True`` when using :meth:`calc_redshift()`.
        """
        return self._z

    @property
    def z_host(self):
        """    
        float or None: 
        The redshift of the localised FRB host galaxy. Note that this an 
        observed quantity, not the estimated redshift :attr:`z` calculated with
        :meth:`calc_redshift()`
        """
        return self._z_host

#    @property
#    def z_uncert(self):
#        """float: The uncertainty in the measurement of the FRB redshift"""
#        return self._z_uncert

#    @property
#    def scatt_index(self):
#        """float: The scattering index of the FRB."""
#        return self._scatt_index

    @property
    def w_obs(self):
        """
        float or None: 
        The observed width of the pulse obtained by a pulse fitting algorithm.
        Units: :math:`\\rm{ms}`
        """
        return self._w_obs

    @property
    def s_peak_obs(self):
        """
        float or None: 
        The observed peak flux density of the burst. Units: 
        """
        return self._s_peak_obs

    @property
    def f_obs(self):
        """The Milky Way component of the dispersion measure."""
        return self._f_obs

    @property
    def raj(self):
        """
        astropy.coordinates.angles.Longitude or None: The right accension in J2000 coordinates of the best
        estimate of the FRB position.
        """
        return self._raj

    @property
    def decj(self):
        """The Milky Way component of the dispersion measure."""
        return self._decj

    @property
    def gl(self):
        """The Milky Way component of the dispersion measure."""
        return self._gl

    @property
    def gb(self):
        """The Milky Way component of the dispersion measure."""
        return self._gb

    @property
    def skycoords(self):
        """
        astropy.coordinates.sky_coordinate.SkyCoord or None: The 
        skycoords of the FRB. This is calculated from either 
        (:attr:`raj`, :attr:`decj`) or (:attr:`gl`, :attr:`gb`).
        """
        return self._skycoords

    @property
    def dm_igm(self):
        """

        """
        return self._dm_igm
