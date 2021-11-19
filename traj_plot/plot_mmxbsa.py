
from data_plot_1D import *


# TODO: make this it's own class, actually should subclass main plot class
# that way same methods can be accessed
########################################################################
##################### 19F MM/PB-GBSA Results Plots #####################
########################################################################
from scipy.stats import ttest_ind
def get_mmxbsa_stats(gsolv, replicates=5):
    """
    From the 5 1µs replicates, calculate the average, std, and ∆∆G (19F-WT).

    Parameters
    ----------
    type : str
        'pb' (default) or 'gb'.
    replicates : int
        Dataset replicates.
    """
    if gsolv == "pb":
        index = 16
        print(f"\nMM/{gsolv.upper()}SA:")
    elif gsolv == "gb":
        index = 15
        print(f"\nMM/{gsolv.upper()}SA:")
    else:
        raise ValueError("gsolv must be 'pb' or 'gb'.")

    systems = ["wt", "w4f", "w5f", "w6f", "w7f"]
    # rows = sys names, cols = average, stdev, ∆G, error
    delta = np.zeros(shape=(len(systems), 4))

    # array of all data, each col = sys, 500 frames of data with 5 replicates
    all_sys_data = np.zeros(shape=(401*5, len(systems)))

    # array of just the averages of each 1µs replicate, each col = sys
    rep_data = np.zeros(shape=(replicates, len(systems)))

    # fill out each row at a time
    for num, sys in enumerate(["wt", "w4f", "w5f", "w6f", "w7f"]):  
        
        data = np.array([pre_processing(data=np.genfromtxt(f"ipq/{sys}/v{i:02d}/1us_noion/MM{gsolv.upper()}SA_ENERGY_OUT_INP1.dat", 
                        skip_header=2, delimiter=","), time_units=500, index=index)[1] for i in range(0, replicates)])
        averages = [np.average(row) for row in data]
        # SEM = std / sqrt(n) : SEM estimates sample mean with greater precsion as n increases
        # standard error of the mean is just the standard deviation of your estimate of the mean
        sems = [np.std(row) / np.sqrt(len(row)) for row in data]

        # there are multiple error propagation steps: SEM of 500 frames of 1us trial, SEM of 5 1us trials, then ∆G error
        if sys == "wt":
            wt_avg = np.average(averages)
            wt_se = np.std(averages) / np.sqrt(len(averages))
            # propagated SE = sqrt(SE_trial1^2 + SE_trial2^2 ...SE_trialN^2) where N = 5
            wt_pse = np.sqrt(np.sum(np.square(sems)))
            
        # fill out columns: avg, stdev, ∆G
        delta[num, 0] = np.average(averages)
        #delta[num, 1] = np.std(averages)

        # SEM = std / sqrt(N)
        #delta[num, 1] = np.std(averages) / np.sqrt(len(averages))

        # pSE of 5 sets of 500 frame average SE
        delta[num, 1] = np.sqrt(np.sum(np.square(sems)))

        delta[num, 2] = np.average(averages) - wt_avg
        # error propagated to ∆G (19F-WT)
        # propagated SE = sqrt(SE_19F^2 + SE_WT^2)
        #delta[num, 3] = np.sqrt((np.std(averages) / np.sqrt(len(averages))**2) + wt_se**2)
        # propagated propagated SE
        delta[num, 3] = np.sqrt((np.sqrt(np.sum(np.square(sems)))**2) + wt_pse**2)

        # build the entire dataset array
        data = np.reshape(data, -1)
        all_sys_data[:,num] = data

        # build only the average of each replicate array
        rep_data[:,num] = averages

    # welches t-test to see if W6F is significantly different from WT
    # welch_all = ttest_ind(all_sys_data[0], all_sys_data[3], equal_var=False)

    # welches t-test to see if the 1µs replicate data is different from W6F to WT
    welch_rep = ttest_ind(rep_data[:,0], rep_data[:,3], equal_var=False)
    print(f"Welches t-test (W6F vs WT), p-value = {welch_rep[1]}")

    df = pd.DataFrame(delta, columns=["Avg G", "±pSE", "∆G(19F-WT)", "±pSE"], index=[i.upper() for i in systems])
    pd.options.display.float_format = '{:,.3f}'.format
    print(df)
        
#get_mmxbsa_stats("pb")

def plot_19F_mmxbsa_dist(type="pb", mechanics="mm", index=16, replicates=5, ax=None, linewidth=3, space=100):
    """
    Plot the PB or GB SA distrubution.
    
    Parameters
    ----------
    type : str
        'pb' (default) or 'gb'.
    mechanics : str
        'mm' (default) or 'qm'. Note QM is for QM/MM GBSA. 
    index : int
        Default 16 for PBSA, which is the total avg G. Other columns are:
        Frame#,BOND,ANGLE,DIHED,UB,IMP,CMAP,VDWAALS,EEL,1-4 VDW,1-4 EEL,EPB,ENPOLAR,EDISPER,G gas,G solv,TOTAL
        For GBSA, use 15 for total avg G. Other columns are:
        Frame#,BOND,ANGLE,DIHED,UB,IMP,CMAP,VDWAALS,EEL,1-4 VDW,1-4 EEL,EGB,ESURF,G gas,G solv,TOTAL
    replicates : int
        Dataset replicates.
    """
    cmap = cm.tab10
    norm = Normalize(vmin=0, vmax=10)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10,7))
    else:
        fig = plt.gca()

    if type == "pb":
        key_index = 2
    elif type == "gb":
        key_index = 1
    # make a dict where key = header_val, val = actual_desc
    terms = np.genfromtxt("MMXBSA_KEY.csv", skip_header=1, delimiter=",", dtype=str)
    term_dict = dict(zip(terms[:,key_index], terms[:,0]))

    for num, sys in enumerate(["wt", "w4f", "w5f", "w6f", "w7f"]):
        if sys == "wt":
            color = "dimgrey"
        else:
            color = cmap(norm(num - 1))
        
        # set up 1D data array of the column of interest
        data = np.array([pre_processing(data=np.genfromtxt(f"ipq/{sys}/v{i:02d}/1us_noion/{mechanics.upper()}{type.upper()}SA_ENERGY_OUT.dat", 
                                 skip_header=2, delimiter=","), time_units=500, index=index)[1] for i in range(0, replicates)])
        data = np.reshape(data, -1)

        # get the title and column names
        with open(f"ipq/wt/v00/1us_noion/{mechanics.upper()}{type.upper()}SA_ENERGY_OUT.dat") as f:
            title = f.readline()
            header = f.readline().strip()
        col_names = header.split(",")

        low_lim = np.min(data) - space
        upp_lim = np.max(data) + space

        #ax.hist(data, color=color, alpha=0.75, bins=50)

        # KDE plot of the 1D total dataset for a single variant
        dist_plot(data, ax=ax, color=color, ylim=(low_lim, upp_lim), linewidth=linewidth)

        # plot formatting
        #ax.set_ylim(0, 0.012)
        ax.set_xlim(low_lim, upp_lim)
        #ax.set_xlim(-5800, -5450)
        ax.set_xlabel(r"$\bar{G}$ (kcal/mol)", labelpad=16)
        #ax.set_xticks(np.arange(-5800, -5400, 50))
        #ax.xaxis.set_minor_locator(MultipleLocator(25))
        ax.tick_params(axis="x", which="major", labelsize=18, pad=4.5)

        # labels
        #ax.set_title(f"{title} {col_names[index]} Contribution")
        ax.set_title("INP=2")
        
        # use labels from header
        #ax.set_xlabel(col_names[index], labelpad=16)

        # use the dict to translate the shorthand terms
        #ax.set_xlabel(term_dict[col_names[index]]) #, fontsize=8, labelpad=8)

        #add_patch(ax, 0.01 + 0.2 * num, -0.435, color, sys.upper() + " CypA", fontsize=17, recspace=-0.01)

        #print(f"\n{sys.upper()}: AVG = {np.mean(data)} ± {np.std(data)}")

    #fig.tight_layout()
    #fig.savefig(f"figures/initial_mm{type}sa.png", dpi=300, transparent=True)

get_mmxbsa_stats("pb")
# index 10 is interesting
#plot_19F_mmxbsa_dist("pb", index=16, replicates=5, space=50)
#plot_19F_mmxbsa_dist("gb", "qm", index=16, replicates=5)
#plt.tight_layout()
#plt.show()
#plt.savefig("figures/total_mmpbsa_inp2.png", dpi=300, transparent=True)

def plot_19F_mmpbsa_decomp_dist(index=7, replicates=5, ax=None, linewidth=3, space=100, res=121):
    """
    Plot the PBSA decomp distrubution.
    
    Parameters
    ----------
    index : int
        PB residue decomp, 8 cols:
        Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL
    replicates : int
        Dataset replicates.
    """
    cmap = cm.tab10
    norm = Normalize(vmin=0, vmax=10)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10,7))
    else:
        fig = plt.gca()

    for num, sys in enumerate(["wt", "w4f", "w5f", "w6f", "w7f"]):
        if sys == "wt":
            color = "dimgrey"
        else:
            color = cmap(norm(num - 1))
        
        # set up 1D data array of the column of interest
        # only 121 is no ending, .dat is 120-122
        # TODO:need to find a better way to seperate total, sidechain, and backbone (each block is 401 frames + 4 header)
        # data = np.array([pre_processing(data=np.genfromtxt( 
        #                 open(f"ipq/{sys}/v{i:02d}/1us_noion/MMPBSA_DECOMP_ENERGY_OUT", "r").readlines()[0 : 4 + (401*1)], 
        #                 skip_header=4, delimiter=","), time_units=500, index=index)[1] for i in range(0, replicates)])
        data = np.array([np.genfromtxt(open(f"ipq/{sys}/v{i:02d}/1us_noion/MMPBSA_DECOMP_ENERGY_OUT.dat", "r").readlines()[0 : 4 + (401*3)], 
                        skip_header=4, delimiter=",") for i in range(0, replicates)])
        # shape is 5, 401*X, 8: get only the residue of interest
        data = data[data[:,:,1] == res]
        #print(np.shape(data)) # now 401*5 rows and 8 cols
        # isolate only 1 residue and change to 1D: second col is residue number
        #data = np.reshape(data[index], -1)

        # get the title and column names TODO:need a better way of getting the 4th line
        with open(f"ipq/w4f/v00/1us_noion/MMPBSA_DECOMP_ENERGY_OUT") as f:
            title = f.readline()
            blanck_line = f.readline()
            decomp_title = f.readline()
            header = f.readline().strip()
        col_names = header.split(",")

        low_lim = np.min(data[:,index]) - space
        upp_lim = np.max(data[:,index]) + space

        # KDE plot of the 1D total dataset for a single variant
        dist_plot(data[:,index], ax=ax, color=color, ylim=(low_lim, upp_lim), linewidth=linewidth)

        # plot formatting
        #ax.set_ylim(0, 0.012)
        ax.set_xlim(low_lim, upp_lim)
        #ax.set_xlabel(r"$\bar{G}$ (kcal/mol)", labelpad=16)
        #ax.set_xticks(np.arange(-5800, -5400, 50))
        #ax.xaxis.set_minor_locator(MultipleLocator(25))
        #ax.tick_params(axis="x", which="major", labelsize=18, pad=4.5)

        # labels
        #ax.set_title(f"{title} {col_names[index]} Contribution")
        
        # use labels from header
        ax.set_xlabel(col_names[index], labelpad=8, fontsize=14)

        #add_patch(ax, 0.01 + 0.2 * num, -0.435, color, sys.upper() + " CypA", fontsize=17, recspace=-0.01)

        #print(f"\n{sys.upper()}: AVG = {np.mean(data)} ± {np.std(data)}")

    #fig.tight_layout()
    #plt.show()
    #fig.savefig(f"figures/initial_mm{type}sa.png", dpi=300, transparent=True)

#plot_19F_mmpbsa_decomp_dist(index=7)
#plt.show()

# with open("ipq/w4f/v00/1us_noion/MMPBSA_DECOMP_ENERGY_OUT", "r") as f:
#     lines = f.readlines()
# # only grab total energy contributions of the 401 frames of data + 4 header lines
# data = np.genfromtxt(lines[:405], skip_header=4, delimiter=",")

def multi_xbsa_plot(type="pb", res=121):
    """
    Plot all contributions, skip index 4,5,6 since 0.
    So PBSA has 16-3=13 plots and GBSA has 15-3=12 plots.
    With PB residue decomp, 8 columns are: 
    Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL
    """
    exclude = [0, 4, 5, 6]
    dims = (3,6)
    if type == "pb":
        max = 17
    elif type == "gb":
        max = 17
    elif type == "decomp":
        max = 8
        exclude = [0, 1, 6]
        dims = (3,3)

    fig, ax = plt.subplots(nrows=dims[0], ncols=dims[1], figsize=(13,7))

    # TODO: it would be more efficient to load the entire array and reference each column
    # currently loading new array for each column

    for num, axis in enumerate(fig.axes):
        if num in exclude or num >= max:
            axis.axis("off")
            continue
        if type == "decomp":
            plot_19F_mmpbsa_decomp_dist(num, 5, axis, 2, space=55, res=res)
        elif type == "pb" or type == "gb":
            plot_19F_mmxbsa_dist(type, "qm", num, 5, axis, 1.5, space=40)

    fig.tight_layout()
    #plt.show()
    fig.savefig(f"figures/all_qmgbsa.png", dpi=300, transparent=True)

#multi_xbsa_plot("decomp", 122)
#multi_xbsa_plot("gb")