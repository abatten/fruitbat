"""
FRB
===
Provides definition of the ``Frb`` class.

"""

from . import estimation


class Frb(object):
    """
    Defines the :class:`~Frb` class in the **fruitbat** package.

    Attributes
    ----------
    name : str
        The name of the Frb object.

    dm : float
        The observed dispersion measure of the FRB without Milky Way
        subtraction.
        Units: :math:`\\rm{pc\ cm^{-3}}`

    dm_uncert : float, optional
        The uncertainty in the dispersion measure.
        Units: :math:`\\rm{pc\\ cm^{-3}}` Default: :math:`0.0`

    dm_galaxy : float, optional
        The modeled contribution to the FRB DM by electrons in the Milky Way.
        Units: :math:`\\rm{pc\\ cm^{-3}}` Default: :math:`0.0`

    dm_excess : float or None, optional
        The DM excess of the FRB over the estimated Galactic DM. If
        `dm_excess` is *None*, then `dm_excess` is calculated by
        :math:`\\rm{DM - DM_{galaxy}}`
        Units: :math:`\\rm{pc\\ cm^{-3}}` Default: *None*

    dm_index : float or None, optional
        The dispersion measure index of the burst :math:`\\alpha` such that
        :math:`\\rm{DM} \\propto \\nu^{-\\alpha}` Default: *None*

    z : float or None, optional
        The redshift of the burst. This assumes that the entire `dm_excess`
        arrives from the IGM and the host galaxy of the FRB and any
        surrounding material contribute nothing to the total DM. This
        should be taken as an upper limit to the bursts true redshift.
        Default: *None*

    z_uncert : float, optional
        The uncertainty in the redshift of the FRB.
        Default: :math:`0.0`

    scatt_index : float or None, optional
        The scattering index (:math:`\\beta`) of the FRB pulse. The 
        scattering index describes how the width (:math:`\\rm{W}`) of the 
        FRB pulse evolves with frequency :math:`\\nu` such that 
        :math:`\\rm{W} \\propto \\nu^{-\\beta}` Default: *None*

    snr : float or None, optional
        The signal-to-noise of the burst.
        Default: *None*

    w_obs : float or None, optional
        The observed width of the pulse obtained by a pulse fitting algorithm.
        Units: ms Default: *None*

    s_peak_obs : float or None, optional
        The observed peak flux density of the burst.
        Units: :math:`\\rm{Jy}` Default: *None*

    f_obs : float or None, optional
        The observed fluence of the FRB. If `f_obs` is *None* and both
        `w_obs` and `s_peak_obs` are not *None* then `f_obs` is
        calculated by :math:`\\rm{W_{obs} \\times S_{peak,obs}}`
        Units: :math:`\\rm{Jy\\ ms}` Default: *None*

    raj : str or None, optional
        The right ascension in J2000 coordinates of the pointing centre
        of the detection beam. This corresponds only to the positioning of
        the beam centre. Default: *None*

    decj : str or None, optional
        The declination in J2000 coordinates of the pointing centre of
        the detection beam. This corresponds only to the positioning of
        the beam centre. Default: *None*
    """

    def __init__(self, name, dm, dm_uncert=0.0, dm_galaxy=0.0, dm_excess=None,
                 dm_index=None, z=None, z_uncert=0.0, scatt_index=None,
                 snr=None, w_obs=None, s_peak_obs=None, f_obs=None, raj=None,
                 decj=None):

        self.name = name
        self.dm = dm
        self.dm_uncert = dm_uncert
        self.dm_galaxy = dm_galaxy

        # Calculate dm_excess from existing parameters if it is not given.
        if not dm_excess:
            self.calc_dm_excess()
        else:
            self.dm_excess = dm_excess

        self.dm_index = dm_index
        self.z = z
        self.z_uncert = z_uncert
        self.scatt_index = scatt_index
        self.snr = snr
        self.w_obs = w_obs
        self.s_peak_obs = s_peak_obs

        # Calculate F_obs if s_peak_obs and w_obs are given
        if (not f_obs) and (s_peak_obs and w_obs):
            self.calc_f_obs()
        else:
            self.f_obs = f_obs
        self.raj = raj
        self.decj = decj

    def __repr__(self):
        return 'Frb({0})'.format(vars(self))

    def calc_redshift(self, method='batten2019', cosmology="planck2018+bao"):
        """
        Calculate the redshift of the FRB from its dm or dm_excess

        Parameters
        ----------
        method : str, optional
            The approximation to use when calculating the redshift. 
            Avaliable methods: ``'batten2019'``, ``zhang2018``,
            ``'inoue2004'``, ``'ioka2003'``

        cosmology : str, optional
            The method `inoue2004` has the option to choose which cosmology 
            to assume when performing the redshift estimation.
            Avaliable cosmologies:

        Returns
        -------
        z : float
            The redshift of the FRB.

        z_uncert : float
            The uncertainty of the redshift estimation.
        """

        z, z_uncert = estimation.calc_redshift(self.dm_excess,
                                               self.dm_uncert, 
                                               method, cosmology)

        self.z = z
        self.z_uncert = z_uncert

        return z, z_uncert

    def calc_dm_excess(self):
        """
        Calculates the dispersion measure excess of the

        .. math ::
            \\rm{DM_{excess} = DM - DM_{galaxy}}

        Returns
        -------
        float
            The dispersion measure excess.
        """
        dm_excess = self.dm - self.dm_galaxy
        self.dm_excess = dm_excess
        return dm_excess

    def calc_f_obs(self):
        """
        Calculates the observed fluence of the FRB. This required `w_obs`
        and `s_peak_obs` to both not be None.

        .. math ::
            \\rm{F_{obs} = W_{obs} \\times S_{peak, obs}}
        """

        if (not self.w_obs) or (not self.s_peak_obs):
            err_msg = ("calc_f_obs requires both w_obs and s_peak_obs "
                       "to not be None")
            raise ValueError(err_msg)

        f_obs = self.w_obs * self.s_peak_obs
        self.f_obs = f_obs
        return f_obs
