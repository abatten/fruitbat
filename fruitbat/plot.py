import matplotlib.pyplot as plt
import numpy as np

from . import estimate
from . import cosmology


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


def create_method_comparison(filename="method_comparison", extension="png",
                             usetex=False):
    """
    Create a plot comparing how estimated redshift changes as a function of
    dispersion measure for each DM-z relation.

    Parameters
    ----------
    filename: string, optional
    The filename of the saved figure. Default: "method_comparison"

    extension: string, optional
        The format to save the figure. e.g "png", "pdf", "esp", etc...
        Default: "png"

    Generates
    ---------
    A figure displaying how estimated redshift changes as a function of
    dispersion measure for each of the different cosmologies.
    """
    set_rc_params(usetex)

    methods = estimate.methods()
    dm_vals = np.linspace(0, 3000, 1000)

    fig = plt.figure(figsize=(8, 8), constrained_layout=True)
    ax = fig.add_subplot(111)

    colours = ["#1b9e77", "#d95f02", "#7570b3"]
    label = [r"$\rm{Ioka 2003}$", r"$\rm{Inoue 2004}$", r"$\rm{Zhang 2018}$"]

    for j, method in enumerate(methods):
        z_vals = np.zeros(len(dm_vals))
        for i, dm in enumerate(dm_vals):
            z_vals[i] = estimate.redshift(dm, method=method)
        ax.plot(dm_vals, z_vals, colours[j], label=label[j], linewidth=4)

    ax.set_xlabel(r"$\rm{DM\ \left[pc \ cm^{-3}\right]}$")
    ax.set_ylabel(r"$\rm{Redshift}$")
    plt.legend(loc='lower right', frameon=False)
    plt.savefig(".".join([filename, extension]))


def create_cosmology_comparison(filename="cosmology_comparison",
                                extension="png", usetex=False):
    """
    Create a plot comparing how the estimated redshift changes as a function
    of dispersion mesure for each cosmology.

    Parameters
    ----------
    filename: string, optional
    The filename of the saved figure. Default: "cosmology_comparison"

    extension: string, optional
        The format to save the figure. e.g "png", "pdf", "esp", etc...
        Default: "png"

    Generates
    ---------
    A figure displaying how estimated redshift changes as a function of
    dispersion measure for each of the different cosmologies.
    """
    set_rc_params(usetex)

    cosmologies = cosmology.builtin().keys()
    dm_vals = np.linspace(0, 3000, 1000)

    fig = plt.figure(figsize=(8, 8), constrained_layout=True)
    ax = fig.add_subplot(111)
    axin = ax.inset_axes([0.05, 0.52, 0.45, 0.45])

    colours = ['#a6cee3', '#1f78b4', '#b2df8a',
               '#33a02c', '#fb9a99', '#e31a1c']
    label = [r"$\rm{WMAP5}$", r"$\rm{WMAP7}$", r"$\rm{WMAP9}$",
             r"$\rm{Planck13}$", r"$\rm{Planck15}$", r"$\rm{Planck18}$"]

    for j, cosmo in enumerate(cosmologies):
        if cosmo == "EAGLE":
            continue
        z_vals = np.zeros(len(dm_vals))
        for i, dm in enumerate(dm_vals):
            z_vals[i] = estimate.redshift(dm, cosmology=cosmo)
        ax.plot(dm_vals, z_vals, colours[j], label=label[j], linewidth=2)
        axin.plot(dm_vals, z_vals, colours[j], linewidth=2)

    ax.set_xlabel(r"$\rm{DM\ \left[pc \ cm^{-3}\right]}$")
    ax.set_ylabel(r"$\rm{Redshift}$")
    plt.legend(loc='lower right', frameon=False)

    axin.set_xlim(2800, 3000)
    axin.set_ylim(3.0, 3.25)
    axin.xaxis.set_tick_params(labelsize=8)
    axin.yaxis.set_tick_params(labelsize=8)

    ax.indicate_inset_zoom(axin)
    plt.savefig(".".join([filename, extension]))
