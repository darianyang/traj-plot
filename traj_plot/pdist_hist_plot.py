"""
Process and plot a heatmap of 2 datasets from standard MD simulations.
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
#sys.path.append("/Users/darian/Google/MBSB/Research/Projects/19F-ff15ipq/umbrella")
sys.path.append("..")
#sys.path.insert(0, "../umbrella/")
from umbrella.plot_phipsi_point_pdb import Get_Dihedrals
#sys.path.append("/Users/darian/Google/MBSB/Research/Projects/19F-ff15ipq/CypA")

# Suppress divide-by-zero in log
np.seterr(divide='ignore', invalid='ignore')

# TODO: clean this up and make it more robust, combine with simple data plot script

#plt.rcParams['figure.figsize']= (10,6)
plt.rcParams.update({'font.size': 14})
plt.rcParams["font.family"]="Sans-serif"
plt.rcParams['font.sans-serif'] = 'Verdana'
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 2.5
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.major.width'] = 2.5

def pre_processing(x_files, y_files, x_loc=1, y_loc=1):
    """
    Parameters
    ----------
    x_files : list of str
        One or multiple dataset paths - 2D array with aux value at each frame.
    y_files : list of str
        One or multiple dataset paths - 2D array with aux value at each frame.
    x_loc : int
        Index of the X ndarray column for the target aux data.
    y_loc : int
        Index of the Y ndarray column for the target aux data.

    Returns
    -------
    X : ndarray
    Y : ndarray
    Z : ndarray
    """
    X = np.concatenate([np.genfromtxt(i)[:,x_loc] for i in x_files])
    Y = np.concatenate([np.genfromtxt(i)[:,y_loc] for i in y_files])
 
    # get rid of nan values: return array without (not) True nan values
    X = X[np.logical_not(np.isnan(X))]
    Y = Y[np.logical_not(np.isnan(Y))]

    # numpy equivalent to: ax.hist2d(c2[:,1], aux)
    hist, x_edges, y_edges = np.histogram2d(X, Y, bins=100)
    # let each row list bins with common y range
    hist = np.transpose(hist)
    # convert histogram counts to pdist (-ln(P(x)))
    hist = -np.log(hist / np.max(hist))
    # get bin midpoints
    midpoints_x = (x_edges[:-1] + x_edges[1:]) / 2
    midpoints_y = (y_edges[:-1] + y_edges[1:]) / 2

    return midpoints_x, midpoints_y, hist

def pdist_plot(X, Y, Z, pmax=5, pmin=0, plot_style="heat", ax=None, cmap="afmhot", savefig=None, **plot_options):
    """
    Parameters
    ----------
    X : ndarray
        1D array of midpoints for histogram data.
    Y : ndarray
        1D array of midpoints for histogram data.
    Z : ndarray
        2D array of the histogram values.
    pmax : int
        Max probability (Z) value.
    pmin : int
        Min probability (Z) value.
    plot_style : str
        'heat', or 'contour'
    ax : mpl axes object
    cmap : str
        mpl colormap string.
    savefig : str
        Path to optionally save the figure output.
    **plot_options : kwargs

    # TODO: add option to plot KDE of hist.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(5,4))
    else:
        fig = plt.gca()

    if plot_style == "heat": # TODO: add contour lines
        plot = ax.pcolormesh(X, Y, Z, cmap=cmap, shading="auto", vmin=pmin, vmax=pmax)
    if plot_style == "contour":
        levels = np.arange(pmin, pmax + 1, 1)
        lines = ax.contour(X, Y, Z, levels=levels, colors="black", linewidths=1)
        plot = ax.contourf(X, Y, Z, levels=levels, cmap=cmap)

   # unpack plot options dictionary
    for key, item in plot_options.items():
        if key == "xlabel":
            ax.set_xlabel(item, weight="bold")
        if key == "ylabel":
            ax.set_ylabel(item, weight="bold")
        if key == "xlim":
            ax.set_xlim(item)
        if key == "ylim":
            ax.set_ylim(item)
        if key == "xticks":
            ax.set_xticks(item)
        if key == "xtick_labels":
            ax.set_xticklabels(item)
        if key == "yticks":
            ax.set_yticks(item)
        if key == "ytick_labels":
            ax.set_yticklabels(item)
        if key == "title":
            ax.set_title(item, fontweight="bold", fontsize=13.5)
        if key == "grid":
            ax.grid(item, alpha=0.5)

    #cbar = fig.colorbar(plot)
    #cbar.set_label(r"$\left[-\ln\,P(x)\right]$")
    #cbar.set_label(r"$\Delta F(\vec{x})\,/\,kT$" + "\n" + r"$\left[-\ln\,P(x)\right]$")
    # TODO: add lines

    #fig.tight_layout()
    if savefig:
        fig.savefig(savefig, dpi=300, transparent=True)
    # else:
    #     plt.show()

    return plot

def multi_cypa_dihedral_plot(type="ramachandran"):
    """
    Plot a multi panel figure of probability distrubutions of CypA trp121 phi/psi. 
    #TODO: add label to bottom instead of title, maybe color code the label, choose a good cmap

    Parameters
    ----------
    type : str
        `ramachandran` (phi, psi) or `janin`(chi1, chi2).
    """
    if type == "ramachandran":
        x = "$\phi$"
        y = "$\psi$"
    elif type == "janin":
        x = "$\chi_{1}$"
        y = "$\chi_{2}$"

    systems = ["wt", "w4f", "w5f", "w6f", "w7f"]
    fig, axes = plt.subplots(nrows=1, ncols=len(systems), figsize=(15, 3.5), sharey=True,
                            gridspec_kw={'width_ratios' : [1,1,1,1,1.25]}
                            )
    for res, ax in zip(systems, axes):
        plot_options = {#"xlabel" : x, 
                        "xlabel" : "$\chi_{1}$", 
                        #"ylabel" : "$\Psi$",
                        "title" : res.upper() + " CypA",
                        "xlim" : (-180, 180),
                        "xticks" : np.arange(-180,270,90),
                        "xtick_labels" : ["", "-90", "", "90", ""],
                        "yticks" : np.arange(-180,270,90),
                        "ytick_labels" : ["-180", "-90", "0", "90", "180"],
                        "ylim" : (-180, 180),
                        "grid" : True,
                        }
        if type == "ramachandran":
            paths = [f"ipq/{res}/v0{i}/1us_noion/dihedral_trp121.dat" for i in range(0,5)]
        elif type == "janin":
            paths = [f"ipq/{res}/v0{i}/1us_noion/dihedral_chi1_2_trp121.dat" for i in range(0,5)]
                
        # phi in second col, psi in third col
        X, Y, Z = pre_processing(paths, paths, x_loc=1, y_loc=2)
        plot = pdist_plot(X, Y, Z, pmax=10, ax=ax, cmap="afmhot_r", **plot_options)

        if res == "wt":
            dh = Get_Dihedrals(f"19F_CYPA_XTALS/{res.upper()}CypA.pdb", "TRP")
        else:
            dh = Get_Dihedrals(f"19F_CYPA_XTALS/{res.upper()}CypA.pdb", res.upper())

        if type == "ramachandran":
            dh.calc_phi_psi()
            print(dh.phi, dh.psi)
            ax.scatter(dh.phi, dh.psi, color="green")
        elif type == "janin":
            dh.calc_chi1_chi2()
            print(dh.chi1, dh.chi2)
            ax.scatter(dh.chi1, dh.chi2, color="green")
        

    cbar = fig.colorbar(plot, pad=0.1)
    cbar.set_label("-ln (P(x))", weight="bold")
    axes[0].set_ylabel(y, weight="bold", labelpad=-1, rotation=0, ha="center", va="center")
    fig.tight_layout()
    plt.show()
    fig.savefig(f"figures/cypa_trp121_{type}_5us.png", dpi=300, transparent=True)

#multi_cypa_dihedral_plot(type="janin")

# plot_options = {"xlabel" : "$\Phi$", 
#                 "ylabel" : "$\Psi$",
#                 "title" : "W4F - WT CypA",
#                 "xlim" : (-180, 180),
#                 "xticks" : np.arange(-180,270,90),
#                 "xtick_labels" : ["", "-90", "", "90", ""],
#                 "yticks" : np.arange(-180,270,90),
#                 "ytick_labels" : ["-180", "-90", "0", "90", "180"],
#                 "ylim" : (-180, 180),
#                 "grid" : True,
#                 }
# # diff plot
# wt_paths = [f"ipq/wt/v0{i}/1us_noion/dihedral_trp121.dat" for i in range(0,5)]
# paths = [f"ipq/w7f/v0{i}/1us_noion/dihedral_trp121.dat" for i in range(0,5)]
# # phi in second col, psi in third col
# X, Y, Za = pre_processing(paths, wt_paths, x_loc=1, y_loc=2)
# X, Y, Zb = pre_processing(paths, paths, x_loc=1, y_loc=2)
# Z = np.subtract(Zb, Za)
# plot = pdist_plot(X, Y, Z, ax=None, cmap="seismic", pmin=-5, pmax=5, **plot_options)
# cbar = plt.colorbar(plot)
# plt.tight_layout()
# plt.show()


def multi_wt_vs_xtal_rms_plot():

    systems = ["w4f", "w5f", "w6f", "w7f"]
    fig, axes = plt.subplots(nrows=1, ncols=len(systems), figsize=(13, 3.5), sharey=True,
                            gridspec_kw={'width_ratios' : [1,1,1,1.25]}
                            )
    for res, ax in zip(systems, axes):
        plot_options = {"xlabel" : "121 RMSD to 3k0n", 
                        #"ylabel" : "121 RMSD to $^{19}F CypA$",
                        "title" : res.upper() + " CypA",
                        "xlim" : (0, 7),
                        #"xticks" : np.arange(2,8,1),
                        # "xtick_labels" : ["", "-90", "", "90", ""],
                        #"yticks" : np.arange(2,8,1),
                        # "ytick_labels" : ["-180", "-90", "0", "90", "180"],
                        "ylim" : (0, 7),
                        "grid" : True,
                        }

        wt_ref_paths = [f"ipq/{res}/v0{i}/1us_noion/rmsd_121_3k0n_ref.dat" for i in range(0,5)]
        xtal_ref_paths = [f"ipq/{res}/v0{i}/1us_noion/rmsd_121_mlwt_ref.dat" for i in range(0,5)]
                
        # phi in second col, psi in third col
        X, Y, Z = pre_processing(wt_ref_paths, xtal_ref_paths, x_loc=1, y_loc=1)
        plot = pdist_plot(X, Y, Z, pmax=10, ax=ax, cmap="afmhot_r", **plot_options)
        # draw diagonal line
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle="--", color="k")

    cbar = fig.colorbar(plot, pad=0.1)
    cbar.set_label("-ln (P(x))", weight="bold")
    axes[0].set_ylabel("121 RMSD to ML non-F", weight="bold")
    fig.tight_layout()
    plt.show()
    fig.savefig(f"figures/cypa_rms_121_3k0n_vs_mlwt_5us.png", dpi=300, transparent=True)

multi_wt_vs_xtal_rms_plot()