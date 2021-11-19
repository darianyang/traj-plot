"""
Calculate the intermolecular native contacts of each 2KOD monomer residue.
Eventually have options to output per-residue time series and pdb file with b-factor or occupancy replaced.

TODO:
Using MDAnalysis, calculate native/non-native contacts.
Map residue pairs involved onto a scatter or heatmap.
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import contacts

def traj_loader(parm, crd, step=1000):
    """
    Load and return a trajectory from mda universe.
    Input parmater file, coordinate file, and loading step interval (default=1000).
    """
    traj = mda.Universe(parm, crd, in_memory=True, in_memory_step=step, verbose=True)
    return traj

def calc_single_residue_contact_residues(universe, selection_1, selection_2, radius=4.5):
    """
    Find the residues between selection_1 and selection_2 within radius.
    TODO: maybe I can use my previous functions and adapt instead?

    Parameters
    ----------
    universe : MDA universe object
        Simulation trajectory.
    selection_1 : str
        Select a single residue.
    selection_2 : str
        Select the entire region to compare to.
    radius : float (optional, default=4.5)
        Contacts within this cutoff (Angstrom) are considered.

    Returns
    -------
    total, native, non-native : list?
        List of residues within radius.

    """

def contacts_within_cutoff(universe, selection_1, selection_2, radius=4.5, plot=False):
    """
    Input universe object, calculate inter-monomer total contacts, plot time series (optional).
    """
    # selections for distance matrix calculation
    mono1 = universe.select_atoms(selection_1)
    mono2 = universe.select_atoms(selection_2)

    timeseries = []
    for ts in universe.trajectory:
        
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(mono1.positions, mono2.positions)

        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])

    df = pd.DataFrame(timeseries, columns=['Frame', '# Contacts'])
    
    if plot is True:
        average_contacts = np.mean(df.iloc[:, 1])
        # plot time series contacts(t)
        fig, ax = plt.subplots()
        ax.plot(df.iloc[:, 0], df.iloc[:, 1])
        ax.set(xlabel='frame', ylabel='number of contacts',
            title='Number of Contacts, average = {:.2f}'.format(average_contacts))

    return df

def fraction_native_contacts(universe, selection_1, selection_2, plot=False, method="soft_cut", ref=None):
    """
    Input universe object, calculate inter-monomer native contacts, plot time series (optional).
    """
    # reference groups: option to use separate pdb file
    if ref is not None:
        mono1 = traj_loader(ref, ref, step=1).select_atoms(selection_1)
        mono2 = traj_loader(ref, ref, step=1).select_atoms(selection_2)
    # reference groups: option to use first frame of trajectory
    elif ref is None:
        mono1 = universe.select_atoms(selection_1)
        mono2 = universe.select_atoms(selection_2)

    # set up analysis of native contacts
    nc = contacts.Contacts(universe, selection=(selection_1, selection_2), refgroup=(mono1, mono2), 
                                method=method, radius=4.5).run()

    # save as 2 column pd df
    nc_df = pd.DataFrame(nc.timeseries, columns=['Frame', 'Fraction Native Contacts'])
    
    #print(nc.contact_matrix)

    if plot is True:
        average_contacts = np.mean(nc.timeseries[:, 1])
        # plot time series q(t)
        fig, ax = plt.subplots()
        ax.plot(nc.timeseries[:, 0], nc.timeseries[:, 1])
        ax.set(xlabel='frame', ylabel='fraction of native contacts', ylim=(0,1),
            title='Native Contacts, average = {:.2f}'.format(average_contacts))
        plt.show()

    return nc_df

def per_residue_contacts_time_series(universe, selection_1, selection_2, datatype, method="soft_cut", ref=None):
    """
    Parameters: universe object, selction terms 1 and 2. Datatype must be 'fraction' (NC) or 'total' (contacts).
    Returns: df of per residue `datatype` contacts time series.
    """
    hxb2_seq = "PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL"
    seq_dict = dict(zip(range(1, 232), [char for char in hxb2_seq]))

    # empty df to be filled in for heatmap
    mono1_df = pd.DataFrame()
    mono2_df = pd.DataFrame()

    # calculate fraction native contacts per residue (monomer 1 - last 20?)
    for res in range(6, 76):
        if datatype == "fraction":
            res_df = fraction_native_contacts(universe, f"resnum {res} and not name H*", selection_2,
                                            method=method, ref=ref)
        elif datatype == "total":
            res_df = contacts_within_cutoff(universe, f"resnum {res} and not name H*", selection_2)

        # canonical numbering correction of + 143 (canonical CTD - 2KOD starts at 144: res 6 = S149)
        mono1_df[seq_dict[res + 143] + str(res + 143)] = res_df.iloc[:, 1] # all rows, second col

    # fill NaN values with 0: no inter monomer native contacts along those residues
    mono1_df = mono1_df.fillna(0).transpose()

    # calculate fraction native contacts per residue (monomer 2 - last 20?)
    for res in range(94, 164):
        if datatype == "fraction":
            res_df = fraction_native_contacts(universe, selection_1, f"resnum {res} and not name H*",
                                            method=method, ref=ref)
        elif datatype == "total":
            res_df = contacts_within_cutoff(universe, selection_1, f"resnum {res} and not name H*")

        # canonical numbering correction: + 143 (canonical CTD) - 88 (monomer 2 numbering)
        mono2_df[seq_dict[res + 143 - 88] + str(res + 143 - 88)] = res_df.iloc[:, 1] # all rows, second col

    # fill NaN values with 0: no inter monomer contacts along those residues
    mono2_df = mono2_df.fillna(0).transpose()

    # returns heatmap ready dataframes for monomer 1 and monomer 2 residues
    return mono1_df, mono2_df

def contacts_heatmap(monomer_df_1, monomer_df_2, datatype):
    """
    Plots a heatmap of both monomer 1 and 2.
    Datatype must be "total" or "fraction" contacts data.
    """
    # format and plot heatmap
    fig, ax = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(12,8),
                            gridspec_kw={'width_ratios' : [20, 20, 1.5]})

    if datatype == "fraction":
        cbar_label = "Fraction of Native Contacts"
    elif datatype == "total":
        cbar_label = "Amount of Total Contacts"
    cbar_color = "viridis"

    a = sns.heatmap(monomer_df_1, ax=ax[0], 
                cmap=cbar_color, cbar=False,
                xticklabels=25, yticklabels=1)
    a.set_yticklabels(a.get_ymajorticklabels(), size=8)
    ax[0].set(xlabel="Frame", ylabel="Residue Index", title="CTD Monomer 1")

    b = sns.heatmap(monomer_df_2, ax=ax[1],
                cmap=cbar_color, cbar_kws={'label': cbar_label}, cbar_ax=ax[2],
                xticklabels=25, yticklabels=1)
    b.set_yticklabels(b.get_ymajorticklabels(), size=8)
    ax[1].set(xlabel="Frame", title="CTD Monomer 2")

    plt.tight_layout()
    plt.show()
    #plt.savefig(f"figures/per_res_{datatype}_contacts.png", dpi=300)

def timeseries_to_average_contacts(dataframe):
    """
    Input: dataframe - rows = residue index, column = time series data.
    Returns: 1 column dataframe rows = residue index, column = average contacts data.
    TODO: prob don't need this as sep function.
    """
    df = pd.DataFrame()
    
    # add new column of mean row (axis=1) values: this seems to perserve row index values
    df["Average Contacts"] = dataframe.mean(axis=1)
    return df

def per_residue_contacts_average(dataframe_m1, dataframe_m2, datatype):
    """
    Make a plot with residue index (x) vs average contact value (plotting function).
    This will be good for comparing to CSP and broadening data.
    Input dataframes for monomer 1 and 2 are formatted: row = residue index, column = time series data.
    """
    df_m1 = timeseries_to_average_contacts(dataframe_m1)
    df_m2 = timeseries_to_average_contacts(dataframe_m2)

    if datatype == "fraction":
        y_label = "Fraction of Native Contacts"
    elif datatype == "total":
        y_label = "Amount of Total Contacts"

    fig, ax = plt.subplots(figsize=(12,8))

    ax.step(df_m1.index, df_m1.iloc[:, 0], where="mid", alpha=0.5, label="CTD Monomer 1")
    ax.step(df_m2.index, df_m2.iloc[:, 0], where="mid", alpha=0.5, label="CTD Monomer 2")

    ax.set_xticks(range(0, len(df_m1.iloc[:, 0]), 1))
    ax.set_xticklabels(df_m1.index, rotation=60, fontsize=8)

    ax.set(ylabel=y_label, xlabel="Residue Index")

    ax.legend(loc=1)
    
    fig.tight_layout()
    plt.show()

# load trajectory
model = "m01"
m01 = traj_loader(f"{model}/{model}_2kod_dry.prmtop", f"{model}/{model}_2kod_200ns_100i_dry.ncdf", step=10)

# selection terms
sel_mono1 = "resnum 6:75 and not name H*"
sel_mono2 = "resnum 94:163 and not name H*"
### include n-termini ###
# sel_mono1 = "resnum 1:75 and not name H*"
# sel_mono2 = "resnum 89:163 and not name H*"

# optional reference structure
min_pdb_file = "/Users/darian/Google/MBSB/Research/Projects/hiv1_capsid_sim/2kod_std_sim/pdb_timepoints/02_min.pdb"
pdb_file = \
"/Users/darian/Google/MBSB/Research/Projects/hiv1_capsid_sim/2kod_std_sim/pdb_timepoints/m01_2kod_amb_main_ref.pdb"


#------------------------------ single residue time series -------------------------------#
# fraction_native_contacts(m01, sel_mono1, sel_mono2, plot=True, method="hard_cut", ref=None)
# contacts_within_cutoff(m01, sel_mono1, sel_mono2, plot=True)
# plt.show()

#--------------------------- all residues time series heatmap ------------------------------#
datatype = "fraction"
df_m1, df_m2 = per_residue_contacts_time_series(m01, sel_mono1, sel_mono2, datatype,
                                                method="soft_cut", ref=min_pdb_file)
#contacts_heatmap(df_m1, df_m2, datatype)

#-------------------------- all residues average contacts (1-D) ----------------------------#
#per_residue_contacts_average(df_m1, df_m2, datatype)


#---------------------------- pdb average contact data output ------------------------------#
from pdb_bfactor_fill import replace_bf_pdb

# original pdb file (2KOD M01)
cannon_pdb_file = "/Users/Darian/Research/Projects/hiv1_capsid_sim/pdb_files/2KOD/m01/m01_2kod_pH_amb.pdb"

# convert time series data to avg contact data, remove column name,
# remove aa 1 letter code from row index, revert canonical naming to linear (1-176)
m1_avg = timeseries_to_average_contacts(df_m1).rename(columns={"Average Contacts":""}, 
                                                      index=lambda x: int(x[1:]) - 143)
m2_avg = timeseries_to_average_contacts(df_m2).rename(columns={"Average Contacts":""}, 
                                                      index=lambda x: int(x[1:]) - 143 + 88)

contacts_out = f"m01/200ns_nc/2kod_{datatype}_contacts"

# combine and write out both monomers to a tsv file
pd.concat([m1_avg, m2_avg]).to_csv(f'{contacts_out}.tsv', sep='\t')

# replace b-factor column of pdb with contact data
replace_bf_pdb(cannon_pdb_file, f'{contacts_out}.pdb', f'{contacts_out}.tsv')

# TODO: could use the min pdb file but need to change the numbering setup in pdb_bfactor_fill
# perhaps add option/arg for naming: linear vs canonical / uniprot style