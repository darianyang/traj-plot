"""
Simple plot of ired data
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize']= (10,6)
plt.rcParams.update({'font.size': 18})
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

def create_relax_arrays(paths):
    """
    From multiple file paths, make R1 and R2 arrays of a single res class.
    """
    r1 = [ 1 / np.genfromtxt(i)[1] for i in paths ]
    r2 = [ 1 / np.genfromtxt(i)[2] for i in paths ]
    return r1, r2

def single_relax_plot(r1, r2, exp, title=None):
    """
    Parameters
    ----------
    r1 : array
        R1 relaxation rates for multiple replicates.
    r2 : array
        R2 relaxation rates for multiple replicates.
    exp : tuple
        Two element tuple of the experimental R1 and R2 rates.

    Returns
    -------
    """
    c1 = "cornflowerblue"
    c2 = "darkorange"
    fig, ax = plt.subplots()
    reps = np.arange(1, len(r1) + 1, 1)
    ax.scatter(reps, r1, label="R1", color=c1)
    ax.scatter(reps, r2, label="R2", color=c2)
    ax.axhline(exp[0], color=c1, linestyle="--")
    ax.axhline(exp[1], color=c2, linestyle="--")
    ax.set_xticks(reps)
    ax.set(xlabel="Replicate", ylabel=r"Relaxation Rate ($s^{-1}$)", title=title, ylim=(0,8.5))
    ax.legend()
    ax.grid(alpha=0.5)
    fig.tight_layout()
    if title:
        fig.savefig(f"figures/{title}_ired.png", dpi=300, transparent=True)
    else:
        fig.savefig("figures/test.png", dpi=300, transparent=False)

def multi_relax_arrays(multi_r1, multi_r2):
    """
    Average and standard deviation arrays of multiple R1 and R2 arrays.
    """
    r1_avg = [np.average(i) for i in multi_r1]
    r1_std = [np.std(i) for i in multi_r1]
    r2_avg = [np.average(i) for i in multi_r2]
    r2_std = [np.std(i) for i in multi_r2]

    return r1_avg, r1_std, r2_avg, r2_std

def multi_relax_plot(r1_avg, r1_std, r2_avg, r2_std, exps, labels):
    """
    Parameters
    ----------

    Returns
    -------
    """
    fig, ax = plt.subplots()
    res_pos = np.arange(1, len(r1_avg) + 1, 1)
    ax.errorbar(res_pos, r1_avg, yerr=r1_std, label="R1", fmt=".", ms=10, elinewidth=0.5)
    ax.errorbar(res_pos, r2_avg, yerr=r2_std, label="R2", fmt=".", ms=10, elinewidth=0.5)

    for i in range(0, len(res_pos)):
        ax.text(res_pos[i], exps[i][0], "-", ha="center", va="center", color="cornflowerblue", fontsize=50, alpha=0.8)
        ax.text(res_pos[i], exps[i][1], "-", ha="center", va="center", color="darkorange", fontsize=50, alpha=0.8)

    ax.set_xticks(res_pos)
    ax.set_xticklabels(labels)
    ax.set(xlabel="Res Class", ylabel=r"Relaxation Rate ($s^{-1}$)", title="All 19F-Trp Data", ylim=(0, 2.5))
    ax.grid(alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/all_ired.png", dpi=300, transparent=False)

# TODO: these are single residue rates... need to update
# exp_rates = {"w4f":(0.99, 1.38), "w5f":(0.75, 0.89), "w6f":(0.67, 0.78), "w7f":(0.87, 1.03)}

# r1_all = []
# r2_all = []
# for key, value in exp_rates.items():
#     r1, r2 = create_relax_arrays([f"ipq/{key}/v{str(i).zfill(2)}/200ns/noe" for i in range(0, 5)])
#     #single_relax_plot(r1, r2, value, title=key)
#     r1_all.append(r1)
#     r2_all.append(r2)

#r1_avg, r1_std, r2_avg, r2_std = multi_relax_arrays(r1_all, r2_all)
#multi_relax_plot(r1_avg, r1_std, r2_avg, r2_std, list(exp_rates.values()), list(exp_rates.keys()))

r1, r2 = create_relax_arrays([f"w4f/v{str(i).zfill(2)}/200ns/noe" for i in range(0, 5)])
single_relax_plot(r1, r2, (0.99, 1.38), title="W4F - G2 REST FRCMOD")