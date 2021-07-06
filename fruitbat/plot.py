import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

from fruitbat import cosmologies, methods, table, utils

def _fruitbat_colors():
    """
    The official list of fruitbat plotting colours.
    """
    color_dict = {
        "black": "#000000",
        "orange": "#FF6800",
        "red": "#C10020",
        "blue": "#007FFF",
    }
    return color_dict

_color_dict = _fruitbat_colors()


def set_rc_params(usetex=False):
    """
    Set the rcParams that will be used in all the plots.
    """

    rc_params = {
        #"axes.prop_cycle": cycler('color',
        #    ['#1b9e77','#d95f02','#7570b3',
        #     '#e7298a','#66a61e','#e6ab02',
        #     '#a6761d','#666666']),
        "axes.labelsize": 18,
        "figure.dpi": 150,
        "legend.fontsize": 12,
        "legend.frameon": False,
        "text.usetex": usetex,
        "xtick.direction": 'in',
        "xtick.labelsize": 14,
        "xtick.minor.visible": True,
        "xtick.top": True,
        "ytick.direction": 'in',
        "ytick.labelsize": 14,
        "ytick.minor.visible": True,
        "ytick.right": True,
    }

    return rc_params





def redshift_pdf(frb, method="Batten2021", sigma=1, usetex=True,
                 filename=None, outputdir=None):
    """
    Plots the redshift pdf and confidence interval for an FRB.

    Parameters
    ----------
    frb : :obj:`fruitbat.Frb`
        An instance of the :obj:`fruitbat.Frb` class.

    method: str, optional

        Default: "Batten2021"

    sigma: [1, 2, 3, 4, 5], optional
        The width of the confidence interval in units of sigma.
        Default: 1

    usetex: bool, optional
        Use LaTeX font when creating the figure. Set this to false to
        disable Latex fonts.
        Default: True

    filename: str or None, optional

    Returns
    -------
    fig: , optional

    ax: , optional

    """
    # Update rcParams for consistent plotting style
    plt.rcParams.update(set_rc_params(usetex=usetex))

    # Calculate the z pdf for the FRB
    zvals, pdf, dz = frb.calc_redshift_pdf(method=method)

    redshift, conf_int = frb.calc_redshift_conf_int(method=method, sigma=sigma)
    conf_int_lower, conf_int_upper = conf_int


    fig, ax = plt.subplots(ncols=1, nrows=1, constrained_layout=True)

    # Plot the redshift pdf
    ax.plot(zvals, pdf, linewidth=2, color=_color_dict["black"])

    # Plot vertical dotted line for the redshift median
    ax.axvline(redshift, linestyle="--", color=_color_dict["blue"],
               linewidth=2)

    conf_int_range = np.linspace(conf_int_lower, conf_int_upper, 100)
    interpolated_z_pdf = interpolate.interp1d(zvals, pdf)

    # Colour in between the confidence interval
    ax.fill_between(conf_int_range, 0, interpolated_z_pdf(conf_int_range),
                    color=_color_dict["orange"], alpha=0.5)

    ax.set_xlim(0, zvals[-1])
    ax.set_ylim(0, 1.05 * np.max(pdf))

    ax.set_xlabel(r"$\mathrm{Redshift}$")
    ax.set_ylabel(r"$P(z | \mathrm{DM}) P(z)$")

    # Set the position of the plot text
    text_ypos_top = 0.90
    text_ypos_padding = 0.06  # The space between each of the lines of text
    if redshift < 1.5:
        text_xpos = 0.60   # Low redshift -> text on right
    else:
        text_xpos = 0.05   # High redshift -> text on left

    text_items = {
        "name"     : (None if frb.name is None
                      else r"$\mathrm{{%s}}$" % frb.name),

        "dm"       : (r"$\mathrm{{DM}} = {}\ \mathrm{{pc\ "
                      r"cm^{{-3}}}}$".format(frb.dm.value)),

        "dm_galaxy": (r"$\mathrm{{DM_{{MW}}}} = {:.1f}\ "
                      r"\mathrm{{pc\ cm^{{-3}}}}$".format(frb.dm_galaxy.value)),

        "dm_model" : (r"$\mathrm{{Galactic\ DM\ Model = %s}}$"
                      % frb.dm_galaxy_model),

        "dm_excess": (r"$\mathrm{{DM_{{Excess}}}} = {:.1f}\ \mathrm{{pc\ "
                      r"cm^{{-3}}}}$".format(frb.dm_excess.value)),

        "redshift" : (r"$z = {:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$".format(redshift,
                      redshift - conf_int_lower, conf_int_upper - redshift)),
    }

    # Some of the text_items may be None, we cant to count the number of
    # non-None items. So I have set up to lune_number to count them.
    # If I don't do this the text is placed in the wrong place when
    # the name is None.
    line_number = 0
    for item in text_items:
        if text_items[item] is not None:
            ax.text(text_xpos,
                    text_ypos_top - text_ypos_padding * line_number, text_items[item],
                    horizontalalignment='left',
                    verticalalignment='center',
                    transform=ax.transAxes,
                    fontsize=11)
            line_number += 1

    if filename is not None:
        if outputdir is not None:
            os.path.join(outputdir, filename)
        plt.savefig("{}.png".format(filename))
        plt.close()

    elif filename is None:
        return fig, ax



def method_comparison(filename=None, extension="png", usetex=False,
                      passed_ax=None, **kwargs):
    """
    Create a plot comparing how estimated redshift changes as a
    function of dispersion measure for each DM-z relation.

    Parameters
    ----------
    filename: string or None, optional
        The filename of the saved figure. Default: *None*

    extension: string, optional
        The format to save the figure. e.g "png", "pdf", "eps", etc...
        Default: "png"

    usetex: bool, optional
        Use LaTeX for for fonts.

    passed_ax: or None, optional

    Generates
    ---------
    A figure displaying how estimated redshift changes as a function of
    dispersion measure for each of the different cosmologies.
    """
    # Update rcParams for consistent plotting style
    plt.rcParams.update(set_rc_params(usetex=usetex))

    if passed_ax:
        ax = passed_ax
    else:
        fig = plt.figure(figsize=(8, 8), constrained_layout=True)
        ax = fig.add_subplot(111)

    method_list = methods.builtin_method_functions()
    method_list.pop("Batten2021")
    dm_vals = np.linspace(0, 3000, 1000)

    colours = ["#1b9e77", "#d95f02", "#7570b3"]
    label = [r"$\rm{Ioka 2003}$", r"$\rm{Inoue 2004}$", r"$\rm{Zhang 2018}$"]

    for j, method in enumerate(method_list):
        z_vals = np.zeros(len(dm_vals))
        if 'cosmology' in kwargs:
            cosmology = kwargs['cosmology']
        else:
            cosmology = 'Planck18'

        table_name = "{}.hdf5".format(method)
        table_name = utils.get_path_to_file_from_here(table_name, subdirs=["data"])


        for i, dm in enumerate(dm_vals):
            z_vals[i] = table.get_z_from_table(dm, table_name, cosmology)

        ax.plot(dm_vals, z_vals, colours[j], label=label[j], **kwargs)

    if not passed_ax:
        ax.set_ylabel(r"$\rm{Redshift}$")


    ax.set_xlabel(r"$\rm{DM\ \left[pc \ cm^{-3}\right]}$")
    ax.legend(loc='lower right', frameon=False)

    if filename is not None:
        plt.savefig(".".join([filename, extension]))

    if passed_ax:
        return ax
    else:
        return fig


def cosmology_comparison(filename="", extension="png", usetex=False,
                         passed_ax=None, **kwargs):
    """
    Create a plot comparing how the estimated redshift changes as a
    function of dispersion mesure for each cosmology.

    Parameters
    ----------
    filename: string, optional
    The filename of the saved figure. Default: "cosmology_comparison"

    extension: string, optional
        The format to save the figure. e.g "png", "pdf", "eps", etc...
        Default: "png"

    Generates
    ---------
    A figure displaying how estimated redshift changes as a function of
    dispersion measure for each of the different cosmologies.
    """
    # Update rcParams for consistent plotting style
    plt.rcParams.update(set_rc_params(usetex=usetex))

    if passed_ax:
        ax = passed_ax
    else:
        fig = plt.figure(figsize=(8, 8), constrained_layout=True)
        ax = fig.add_subplot(111)

    # Remove EAGLE from cosmologies since it is the same as Planck13
    cosmology_list = cosmologies.builtin_cosmology_functions()
    cosmology_list.pop("EAGLE")

    dm_vals = np.linspace(0, 3000, 1000)

    add_axin = True
    try:
        # Add inset plot showing the part where cosmologies diverge the most.
        axin = ax.inset_axes([0.05, 0.52, 0.45, 0.45])
    except Exception:
        print("""Skipping inset axis in cosmology plot. Requires Python 3 and
              Matplotlib 3.0""")
        add_axin = False

    colours = ['#a6cee3', '#1f78b4', '#b2df8a',
               '#33a02c', '#fb9a99', '#e31a1c']
    label = [r"$\rm{WMAP5}$", r"$\rm{WMAP7}$", r"$\rm{WMAP9}$",
             r"$\rm{Planck13}$", r"$\rm{Planck15}$", r"$\rm{Planck18}$"]

    for j, cosmo in enumerate(cosmology_list):
        z_vals = np.zeros(len(dm_vals))
        if 'method' in kwargs:
            method = kwargs['method']
        else:
            method = 'Inoue2004'

        table_name = "{}.hdf5".format(method)
        table_name = utils.get_path_to_file_from_here(table_name, subdirs=["data"])
        for i, dm in enumerate(dm_vals):
            z_vals[i] = table.get_z_from_table(dm, table_name, cosmo)

        ax.plot(dm_vals, z_vals, colours[j], label=label[j], **kwargs)
        if add_axin:
            axin.plot(dm_vals, z_vals, colours[j], **kwargs)

    ax.set_xlabel(r"$\rm{DM\ \left[pc \ cm^{-3}\right]}$")

    if not passed_ax:
        ax.set_ylabel(r"$\rm{Redshift}$")

    ax.legend(loc='lower right', frameon=False)

    if add_axin:
        axin.set_xlim(2800, 3000)
        axin.set_ylim(3.0, 3.25)
        axin.xaxis.set_tick_params(labelsize=8)
        axin.yaxis.set_tick_params(labelsize=8)

        ax.indicate_inset_zoom(axin)

    if filename != "":
        plt.savefig(".".join([filename, extension]))

    if passed_ax:
        return ax
    else:
        return fig
