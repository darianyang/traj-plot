"""
Use cpptraj contacts output to plot native contact information with:
1) heat map of residue contact pair fractions
2) time series of % native contacts over simulation time
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def residue_contact_map(path, cbar_title, diff=None, heatmap=False):
    """
    Build residue contact map from cpptraj res pair output.
    Additional option to generate difference plot of 2 res pair outputs.
    """
    res_pair_frac = np.genfromtxt(path) # also can use np.loadtext()
    # x = monomer 1, y = monomer 2, z = fraction native contacts
    x = res_pair_frac[:, :1]
    y = np.subtract(res_pair_frac[:, 1:2], 88)
    #y = res_pair_frac[:, 1:2]
    z = res_pair_frac[:, 2:3]

    # option for diff contact map
    if diff != None:
        z1 = np.genfromtxt(diff)[:, 2:3]
        z = np.subtract(z, z1)
        plt.title('First 50ns - Last 50ns Dimer Int')

    # heatmap option with matrix input (TODO: dosen't really work)
    if heatmap is True:
        # build pandas df
        df = pd.DataFrame(res_pair_frac)
        df.columns = ['res1', 'res2', 'Strength NC'] #optional step
        print(df)

        shape = int(np.sqrt(len(res_pair_frac)))

        # make matrix of size amount of residue comparisons: CTD dimer = 176 residues
        contact_matrix = np.asarray(df['Strength NC']).reshape(shape, shape)
        
        hm = plt.imshow(contact_matrix, cmap='afmhot_r')
    else:
        # scatter plot with squares and colorbar free energy
        hm = plt.scatter(x, y, c=z, cmap='viridis', marker='s', s=15)

    if diff != None:
        plt.colorbar(hm).set_label('Fraction of total native contacts Difference')
        plt.clim(-5,5)
    else:
        plt.colorbar(hm).set_label(cbar_title)

    plt.xlabel('CTD Monomer 1')
    plt.ylabel('CTD Monomer 2')

    # helix 3
    plt.axvspan(35, 50, color="red", alpha=0.10, label="Helix 3")
    plt.axhspan(35, 50, color="red", alpha=0.10)

    # helix 1
    plt.axvspan(6, 10, color="blue", alpha=0.10, label="Helix 1")
    plt.axhspan(6, 10, color="blue", alpha=0.10)

    #helix_elements = {6:10, 17:30, 35:50, 54:63, 67:75}

    # highlight secondary elements
    # for key, value in helix_elements.items():
    #     plt.axvspan(key, value, color="red", alpha=0.10)
    #     plt.axvspan(key + 88, value + 88, color="red", alpha=0.10)
    #     plt.axhspan(key, value, color="red", alpha=0.10)
    #     plt.axhspan(key + 88, value + 88, color="red", alpha=0.10)

    #plt.clim(0,30)
    plt.legend(loc=1)

    #plt.show()
    #plt.savefig('figures/1us_nc_intra.png', dpi=300)


def percent_native_time_series(path, color, label, total_contacts, nnc=False):
    """
    Plot time series of fraction nc or nnc over time.
    """
    contacts_file = np.genfromtxt(path)
    contacts_file = contacts_file[:200:1] # edit this interval for more sparse dataset
    frame = np.divide(contacts_file[:, :1], 100) # set to proper x-axis unit timestep
    #frame = contacts_file[:, :1]
    contacts = contacts_file[:, 1:2]

    # option to overwrite nc with nnc
    if nnc is True:
        contacts = contacts_file[:, 2:3]

    # fraction of the total amount of contacts, nc or nnc
    percent_contacts = np.divide(contacts, total_contacts)

    plt.plot(frame, percent_contacts, linewidth=1, label=label, color=color)

    plt.xlabel('Time (ns)')
    plt.ylabel('Fraction of Contacts')
    plt.ylim(0, 1)

    #plt.legend(loc=1)

def get_total_contacts(res_pair_file, with_nnc=False):
    """
    Take res pair file and extract the amount of native contacts for the dataset.
    Also option to use 'with nnc=True' which allows extraction of total nc and nnc values.
    Returns int: amount of nc, optionally: amount of nnc
    """
    total_contacts = open(res_pair_file, "r")
    # skip first line (header)
    total_contacts = total_contacts.readlines()[1:]

    # list of native contact amounts only (column 4), append first part of file
    total_nc = []
    for line in total_contacts:
        if '#' not in line:
            total_nc.append(int(line.split()[3]))
        else:
            break
    #print(np.sum(total_nc))

    # additional list of non-native contact amounts, append second part of file
    if with_nnc is True:
        total_nnc = []
        read_nnc = False
        for line in total_contacts:
            if read_nnc is True:
                total_nnc.append(int(line.split()[3])) 
            if '#' in line:
                read_nnc = True

        return np.sum(total_nc), np.sum(total_nnc)
    else:
        return np.sum(total_nc)


#residue_contact_map('m01/1us_nc/native.m01_2kod_nc_map_1us_intra.dat', 'Fraction of total native contacts', heatmap=True)
#residue_contact_map('m05/c3/native.m05_2kod_c3_map_inter.dat', 
#    'Fraction of total native contacts', heatmap=True)

#residue_contact_map('m01/1us_nc/m01_2kod_nc_res_pairs_f50.dat', diff='m01/1us_nc/m01_2kod_nc_res_pairs_l50.dat')
#residue_contact_map('m01/1us_nc/m01_2kod_nc_res_pairs_1us_inter.dat')


# ------------------------------ eq and prod 6.0A time series plot ---------------------------------- #
# nc_inter, nnc_inter = get_total_contacts('m01/equil_1us_nc/m01_2kod_nc_res_pairs_eq_1us_inter_min.dat', with_nnc=True)
# nc_intra, nnc_intra = get_total_contacts('m01/equil_1us_nc/m01_2kod_nc_res_pairs_eq_1us_intra_min.dat', with_nnc=True)

# plt.figure(figsize=(8,6))
# percent_native_time_series('m01/equil_1us_nc/m01_2kod_nc_number_eq_1us_inter_min.dat', 
#     'cornflowerblue', 'Inter Monomer NC', total_contacts=nc_inter)
# percent_native_time_series('m01/equil_1us_nc/m01_2kod_nc_number_eq_1us_inter_min.dat', 
#     'blue', 'Inter Monomer NNC', total_contacts=nnc_inter, nnc=True)
# percent_native_time_series('m01/equil_1us_nc/m01_2kod_nc_number_eq_1us_intra_min.dat', 
#     'darkorange', 'Intra Monomer NC', total_contacts=nc_intra)
# percent_native_time_series('m01/equil_1us_nc/m01_2kod_nc_number_eq_1us_intra_min.dat', 
#     'red', 'Intra Monomer NNC', total_contacts=nnc_intra, nnc=True)
# plt.show()
#plt.savefig('figures/intra_inter_equil_1us_time_series_10i.png', dpi=300)


# ------------------------------ 1us prod 4.5A time series/contact plots ---------------------------------- #
# nc_inter, nnc_inter = get_total_contacts('m01/1us_nc_4.5A/m01_2kod_nc_res_pairs_1us_inter_4.5A.dat', with_nnc=True)
# nc_intra, nnc_intra = get_total_contacts('m01/1us_nc_4.5A/m01_2kod_nc_res_pairs_1us_intra_4.5A.dat', with_nnc=True)

residue_contact_map('m01/1us_nc_4.5A/m01_2kod_nnc_res_pairs_1us_inter_4.5A.dat', 
                    'Amount of total contacts', heatmap=False)
#residue_contact_map('m01/1us_nc_4.5A/m01_2kod_nc_res_pairs_1us_intra_4.5A.dat', 
#                    'Amount of total contacts', heatmap=True)

#plt.show()
plt.savefig("figures/4.5A_1us_nnc_res_pairs.png", dpi=300)


# plt.figure(figsize=(8,6))
# percent_native_time_series('m01/1us_nc_4.5A/m01_2kod_nc_number_1us_inter_4.5A.dat', 
#     'cornflowerblue', 'Inter Monomer NC', total_contacts=nc_inter)
# percent_native_time_series('m01/1us_nc_4.5A/m01_2kod_nc_number_1us_inter_4.5A.dat', 
#     'blue', 'Inter Monomer NNC', total_contacts=nnc_inter, nnc=True)
# percent_native_time_series('m01/1us_nc_4.5A/m01_2kod_nc_number_1us_intra_4.5A.dat', 
#     'darkorange', 'Intra Monomer NC', total_contacts=nc_intra)
# percent_native_time_series('m01/1us_nc_4.5A/m01_2kod_nc_number_1us_intra_4.5A.dat', 
#     'red', 'Intra Monomer NNC', total_contacts=nnc_intra, nnc=True)
# plt.show()
# #plt.savefig('figures/intra_inter_4.5A_time_series_frame1_ref.png', dpi=300)


# ------------------------------ 10ns prod 4.5A time series plot ---------------------------------- #
nc_inter, nnc_inter = get_total_contacts('m01/1us_nc_4.5A/m01_2kod_nc_res_pairs_1us_inter_4.5A.dat', 
                                        with_nnc=True)

# plt.figure(figsize=(8,6))
# percent_native_time_series('m01/1us_nc_4.5A/m01_2kod_nc_number_10ns_inter_4.5A.dat', 
#     'cornflowerblue', 'Inter Monomer NC', total_contacts=nc_inter)

# #plot vertical line for each frame/timepoint that I am observing
# #plt.vlines([0, 0.06, 0.94], 0, 1, colors=['tan', 'turquoise', 'plum'], 
# #           linestyles='dashed', alpha=0.8, label='Timepoints of Interest')

# plt.legend(loc=1)
# plt.show()
#plt.savefig('figures/inter-2ns-timepoints.png', dpi=300)