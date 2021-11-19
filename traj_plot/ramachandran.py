"""
Process and plot pmf data file from wham-2d calculation.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.patheffects as PathEffects
from scipy.interpolate import griddata

#plt.rcParams['figure.figsize']= (8,5)
plt.rcParams.update({'font.size': 12})
plt.rcParams["font.family"]="Sans-serif"
plt.rcParams['font.sans-serif'] = 'Veranda'
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 2.5
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.major.width'] = 2.5

def process_wham(file, bins=72):
    """
    Parameters
    ----------
    file : str
        Path to PMF data file from WHAM output.
    bins : int
        The amount of x and y bins, z will be equal to an array of (bins * bins).

    Returns
    -------
    XY : ndarray
        1 dimensional array of the dihedral range partitioned into n bins. 
    Z : ndarray
        2 dimensional array of (bins * bins) size containing the PMF values at
        each cooresponding dihedral value.
    """
    # extract PMF values
    data = np.genfromtxt(file)
    
    # generate XY array from the PMF binning scheme
    XY = data[0:bins, 1]

    # loop through data array to build 2D Z array
    Z = np.zeros(shape=(bins, bins))
    for dataset in data:
        # populate with PMF at the XY array index that matches the cooresponding X and Y value
        Z[np.where(XY == dataset[0]), np.where(XY == dataset[1])] = dataset[2]
    
    # transpose Z to correct orientation
    return XY, Z.T

def plot_ramachandran(XY, Z, Zdiff=None, ax=None, pmax=5, cmap="afmhot", 
                      title=None, cbar=False, state_labels=False):
    """ 
    TODO: add plot kwargs like with h5_plot?

    Parameters
    ----------
    XY : ndarray
        1D binning scheme for the x and y axes.
    Z : ndarray
        2D heatmap array of the PMF.
    Zdiff : ndarray
        Optional comparison 2D heatmap array to subtract from PMF (Z).
    ax : mpl axis object
        In the figure, axes to plot on. If None, use current axes.
    pmax : int
        Max probability/PMF value, used to set max contour and cbar limits.
        Units are ∆G (kcal/mol).
    cmap : str
        Colormap option.
    title : str
        Optional plot title.
    cbar : bool
        Optionally plot labeled colorbar.
    state_labels : bool
        Optionally plot state labels.

    Returns
    -------
    """
    if Zdiff is not None:
        # hide distracting large differences in depopulated regions
        Z[Z > pmax] = np.nan
        # make difference plot data
        Z = np.subtract(Z, Zdiff)
        # set new lower range
        pmin = -5
    else:
        pmin = 0
    if ax is None:
        fig, ax = plt.subplots(figsize=(4,3.5))
    else:
        fig = plt.gcf()

    levels = np.arange(pmin, pmax + 1, 1)

    # contour or heatmap TODO
    #plot = ax.contourf(XY, XY, Z, levels=levels, cmap=cmap, zorder=0)
    plot = ax.pcolormesh(XY, XY, Z, norm=Normalize(vmin=pmin, vmax=pmax), cmap=cmap, shading="auto", zorder=0)
    lines = ax.contour(XY, XY, Z, levels=levels, colors="black", linewidths=1, zorder=1)
    ax.grid(alpha=0.5, zorder=2)
    
    ticks = (-180, -90, 0, 90, 180)
    tick_labels = [" ", "-90", " ", "-90", " "]

    ax.set_xlabel(r"$\Phi$", labelpad=4, fontweight="bold", fontsize=16)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_ylabel(r"$\Psi$", rotation=1, va='center', ha='left', 
                  labelpad=12, fontweight="bold", fontsize=16)
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels)
    if title:
        ax.set_title(title, fontweight="bold")

    if cbar is True:
        cbar = fig.colorbar(plot)
        cbar.set_label("$\Delta$G (kcal mol$^{-1}$)", fontweight="bold", labelpad=20, rotation=90)
        cbar.add_lines(lines)
        cbar.lines[-1].set_linewidth(2.5)
        cbar.ax.tick_params(size=0)

    states = [("β",       -151,  151),
              ("PPII",     -66,  140),
              ("ξ",       -145,   55),
              ("γ'",       -81,   65),
              ("α",        -70,  -25),
              ("$L_α$",     55,   45),
              ("γ",         73,  -35),
              ("PPII'",     56, -124),
              ("plateau", -100, -130)]
    if state_labels is True:
        for state in states:
            label = ax.text(state[1], state[2], state[0], ha="center", va="center", weight="bold", size=10)
            label.set_path_effects([PathEffects.withStroke(linewidth=1, foreground="w")])

    fig.tight_layout()
    # TODO: savefig option

    return plot, lines

# XY, Z = process_wham(f"PMFs/TRP_pmf.dat")
# plot, lines = plot_ramachandran(XY, Z, title="Legend", state_labels=True, cbar=True)
# plt.savefig("figures/legend_trp.png", dpi=300, transparent=True)
#plt.show()

def format_cbar(grid_axis, plot, lines, label):
    cbar = grid_axis.cax.colorbar(plot)
    cbar.add_lines(lines)
    cbar.lines[-1].set_linewidth(2.5)
    cbar.ax.tick_params(size=0)
    cbar.set_label(label, weight="bold", size=14, labelpad=10)
    for t in cbar.ax.get_yticklabels():
        t.set_horizontalalignment("right")
        t.set_x(3)

def multi_pmf_plot():
    """
    Mutliple panel plot.
    """
    systems = ["TRP", "W4F", "W5F", "W6F", "W7F"]
    #systems = ["TRP", "TYR", "PHE"]
    #systems = ["TYR", "Y3F", "YDF"]
    #systems = ["PHE", "F4F", "FTF"]
    #systems = ["W4F", "W5F", "W6F", "W7F", "Y3F", "YDF", "F4F", "FTF"]
    #systems = ["Y3F", "YDF", "F4F", "FTF"]

    # Set up figure and image grid
    #fig = plt.figure(figsize=(10, 6))
    fig = plt.figure(figsize=(13, 6))

    # grid is single mpl subplot (111) with multiple objects in itself
    grid = ImageGrid(fig, 111,
                     nrows_ncols=(2,len(systems)),
                     axes_pad=(0.15, 0.5),
                     share_all=True,
                     cbar_location="right",
                     cbar_mode="edge",
                     cbar_size="5%",
                     cbar_pad=0.20,
                     )

    # add phi/psi data to image grid
    for ax, res in zip(grid[:len(systems)], systems):
        XY, Z = process_wham(f"PMFs/{res}_pmf.dat")
        plot, lines = plot_ramachandran(XY, Z, ax=ax, pmax=5, title=res, state_labels=False)
        if res == "TRP" or res == "TYR" or res == "PHE":
            continue
        else:
            pdb_data = np.genfromtxt(f"19F_pdbs/dh_data/{res.lower()}_phipsi.tsv", skip_header=1)
            ax.scatter(pdb_data[:, 0], pdb_data[:, 1], color="green", s=10)
        
    # add phi/psi diff data to image grid
    for ax, res in zip(grid[len(systems):], systems):
        XY, Z = process_wham(f"PMFs/{res}_pmf.dat")
        XY, Zdiff = process_wham(f"PMFs/{systems[0]}_pmf.dat")
        plot_diff, lines_diff = plot_ramachandran(XY, Z, Zdiff, ax=ax, pmax=5, 
                                title=f"{res} - {systems[0]}", cmap="seismic")

    format_cbar(grid[0], plot, lines, "$\Delta$G (kcal mol$^{-1}$)")
    format_cbar(grid[len(systems)], plot_diff, lines_diff, "$\Delta$$\Delta$G (kcal mol$^{-1}$)")

    fig.tight_layout()
    #plt.show()
    fig.savefig("figures/trp_diff_bold.png", dpi=300, transparent=True)

multi_pmf_plot()

