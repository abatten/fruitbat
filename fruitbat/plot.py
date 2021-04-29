from __future__ import print_function, absolute_import, division

import matplotlib.pyplot as plt
import numpy as np

from fruitbat import cosmologies, methods, table


def set_rc_params(usetex=False):
    """
    Set the rcParams that will be used in all the plots.
    """
    plt.rcParams["text.usetex"] = usetex
    plt.rcParams["axes.labelsize"] = 20
    plt.rcParams["axes.labelpad"] = 3.0
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["xtick.labelsize"] = 16
    plt.rcParams["ytick.labelsize"] = 16
    plt.rcParams["legend.fontsize"] = 16


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
    set_rc_params(usetex)

    if passed_ax:
        ax = passed_ax
    else:
        fig = plt.figure(figsize=(8, 8), constrained_layout=True)
        ax = fig.add_subplot(111)

    method_list = methods.available_methods()
    dm_vals = np.linspace(0, 3000, 1000)

    colours = ["#1b9e77", "#d95f02", "#7570b3"]
    label = [r"$\rm{Ioka 2003}$", r"$\rm{Inoue 2004}$", r"$\rm{Zhang 2018}$"]

    for j, method in enumerate(method_list):
        z_vals = np.zeros(len(dm_vals))
        if 'cosmology' in kwargs:
            cosmology = kwargs['cosmology']
        else:
            cosmology = 'Planck18'

        table_name = "".join(["_".join([method, cosmology]), ".npz"])
        lookup_table = table.load(table_name)

        for i, dm in enumerate(dm_vals):
            z_vals[i] = table.get_z_from_table(dm, lookup_table)

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
    set_rc_params(usetex)

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

        table_name = "".join(["_".join([method, cosmo]), ".npz"])
        lookup_table = table.load(table_name)        
        for i, dm in enumerate(dm_vals):
            z_vals[i] = table.get_z_from_table(dm, lookup_table)

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

    if filename is not "":
        plt.savefig(".".join([filename, extension]))

    if passed_ax:
        return ax
    else:
        return fig
