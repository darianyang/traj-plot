"""
Take multiple PDB files and build then plot the RMSD matrix.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import MDAnalysis as mda

from pathlib import Path
from MDAnalysis.analysis import rms
from matplotlib import pyplot as plt

# make list of target pdb paths, then load these on the fly for rmsd calc
# 2KOD m01-30
nmr_pdb_paths = []
for pdb in range(1,10):
   nmr_pdb_paths.append(f"/Users/darian/Google/MBSB/Research/Projects/hiv1_capsid_sim/pdb_files/2KOD/m0{pdb}/m0{pdb}_2kod_pH_amb.pdb")
for pdb in range(10,31):
   nmr_pdb_paths.append(f"/Users/darian/Google/MBSB/Research/Projects/hiv1_capsid_sim/pdb_files/2KOD/m{pdb}/m{pdb}_2kod_pH_amb.pdb")

main_pdb_paths = []
main_pdb_paths.append(f"pdb_timepoints_and_ref/m01_2kod_leap.pdb")
#main_pdb_paths.append(f"pdb_timepoints_and_ref/2kod_min.pdb")
main_pdb_paths.append(f"pdb_timepoints_and_ref/1a43_leap.pdb")
#main_pdb_paths.append(f"pdb_timepoints_and_ref/1a43_min.pdb")

# km4 pdb files for m01-05 representative crd
km4_paths = [f"m0{pdb}/m0{pdb}_rep.c{rep}.pdb" 
            for pdb in range(1,6) #m01-05
                for rep in range(0,4) #km=4
            ]

# monomer1 km pdb files for m01-05 representative crd
mono1_paths = [f"m0{pdb}/monomer1/m0{pdb}_ctd1_rep.c{rep}.pdb" 
            for pdb in range(1,6) #m01-05
                for rep in range(0,4) #km=4
            ]

# monomer2 km pdb files for m01-05 representative crd
mono2_paths = [f"m0{pdb}/monomer2/m0{pdb}_ctd2_rep.c{rep}.pdb" 
            for pdb in range(1,6) #m01-05
                for rep in range(0,4) #km=4
            ]

# M2-E175_M1_W184 Dist KM3 2KOD M01 1us
m2_m1_m01_paths = [f"m01/km_dist/m2-m1_rep.c{clust}.pdb"
                   for clust in range(0,3) #km=3
                   ]

# C2 Angle KM3 2KOD M01 1us C2_VEC: 46-49 37-40  
c2_m01_paths_v0 = [f"m01/km_c2_46-49_37-40/m2-m1_rep.c{clust}.pdb"
                   for clust in range(0,3) #km=3
                   ]

# C2 Angle KM3 2KOD M01 1us C2_VEC: 41-44 38-41  
c2_m01_paths_v1 = [f"m01/km_c2_41-44_38-41/m2-m1_rep.c{clust}.pdb"
                   for clust in range(0,3) #km=3
                   ]

# C2 Angle KM3 2KOD M01 1us C2_VEC: 1-75 39
c2_m01_paths_v2 = [f"m01/km_c2_1-75_39/m2-m1_rep.c{clust}.pdb"
                   for clust in range(0,3) #km=3
                   ]     

# C2 Angle KM3 2KOD WE i150 C2_VEC: 1-75 39
c2_we_2kod_paths = [f"agg_2kod_we_v02/i150/km_c2_1-75_39/m2-m1_rep.c{clust}.pdb"
                   for clust in range(0,3) #km=3
                   ]     

# C2 Angle KM3 1A43 WE i150 C2_VEC: 1-75 39
c2_we_1a43_paths = [f"agg_1a43_we_v01/i150/km_c2_1-75_39/m2-m1_rep.c{clust}.pdb"
                   for clust in range(0,3) #km=3
                   ]        

# KM3 2KOD WE i150 M2-M1 Distance
m2_m1_we_2kod_paths = [f"agg_2kod_we_v02/i150/km_dist_m2-m1/m2-m1_rep.c{clust}.pdb"
                      for clust in range(0,3) #km=3
                      ]    

# C2 Angle KM3 1A43 WE i150 M2-M1 Distance
m2_m1_we_1a43_paths = [f"agg_1a43_we_v01/i150/km_dist_m2-m1/m2-m1_rep.c{clust}.pdb"
                      for clust in range(0,3) #km=3
                      ]  

def rmsd_diff_calc(pdb, ref_pdb):
    """
    Take two pdb file path strings, allign and return rmsd value in Angstroms.
    """
    protein = mda.Universe(Path(pdb))
    protein_ref = mda.Universe(Path(ref_pdb))

    # calc non-termini heavy atom rmsd of 2kod
    rmsd = rms.rmsd(#protein.select_atoms('resnum 6:75 and name CA, C, O, N', 'resnum 94:163 and name CA, C, O, N').positions,
                    protein.select_atoms('resnum 6:75 and not name H*', 'resnum 94:163 and not name H*').positions,
                    #protein.select_atoms('resnum 149:218 and name CA, C, O, N').positions, 
                    #protein.select_atoms('resnum 149:218 and not name H*').positions, 
                    #protein_ref.select_atoms('resnum 6:75 and name CA, C, O, N', 'resnum 94:163 and name CA, C, O, N').positions,
                    protein_ref.select_atoms('resnum 6:75 and not name H*', 'resnum 94:163 and not name H*').positions,
                    #protein_ref.select_atoms('resnum 149:218 and not name H*').positions,  
                    center=True,
                    superposition=True)
    return rmsd

def build_pdb_rmsd_matrix(pdb_paths, pdb_diff_path=None):
    """
    Returns rmsd difference matrix for multiple pdb files.
    Returns rmsd_list (3-item list), pdb_comp_amount (int).
    Optional with pdb_diff_path return pdb_diff_comp(int).
    """
    # make 3 column list or ndarray for x, y = (pdb1-n * pdb1-n) and z = rmsd diff
    rmsd_list = [[], [], []]

    # get rmsd difference between each pdb file in nested loop and append
    for pdb0 in pdb_paths:

        # compare 2 different sets of pdb files
        if pdb_diff_path != None:
            for pdb1 in pdb_diff_path:
                # append to x (col 0) pdb in outer loop
                rmsd_list[0].append(pdb_paths.index(pdb0) + 1)
                # append to y (col 1) pdb in inner loop
                rmsd_list[1].append(pdb_diff_path.index(pdb1) + 1)
                
                # find and append to z (col 2) rmsd value between pdb0 and pdb1
                rmsd = rmsd_diff_calc(pdb0, pdb1)
                #print(f"\n    For PDB-A = {pdb0} and PDB-B = {pdb1} : RMSD = {rmsd}")
                rmsd_list[2].append(rmsd)

        elif pdb_diff_path == None:
            for pdb1 in pdb_paths:
                # append to x (col 0) pdb in outer loop
                rmsd_list[0].append(pdb_paths.index(pdb0) + 1)
                # append to y (col 1) pdb in inner loop
                rmsd_list[1].append(pdb_paths.index(pdb1) + 1)
                
                # find and append to z (col 2) rmsd value between pdb0 and pdb1
                rmsd = rmsd_diff_calc(pdb0, pdb1) 
                rmsd_list[2].append(rmsd)

    # amount of pdb files to compare to each other
    pdb_comp_amount = len(pdb_paths)

    if pdb_diff_path == None:
        return rmsd_list, pdb_comp_amount
    elif pdb_diff_path !=None:
        pdb_diff_comp = len(pdb_diff_path)
        return rmsd_list, pdb_comp_amount, pdb_diff_comp

def plot_pdb_rmsd_matrix(rmsd_list, pdb_amount, pdb_diff_amount=None, xylabels=None):
    """
    Generate and plot heatmap of rmsd matrix from list: x*y=pdb z=rmsd. 
    """
    # need to properly format data and build heatmap using numpy
    # tutorials online: e.g. https://blog.quantinsti.com/creating-heatmap-using-python-seaborn/

    x_name = 'Structure' 
    y_name = 'Structures'

    # build pandas df
    df = pd.DataFrame(rmsd_list).transpose()
    df.columns = [x_name, y_name, 'rmsd'] #optional step
    #print(df)

    # TODO: perhaps some adjustment of row and col labels to better match

    # made matrix size dynamic as the len of pdb_path list
    # TODO: I can just use int(np.sqrt(len(df))) as shape
    if pdb_diff_amount == None:
        rmsd_matrix = np.asarray(df['rmsd']).reshape(pdb_amount, pdb_amount)
    elif pdb_diff_amount != None:
        rmsd_matrix = np.asarray(df['rmsd']).reshape(pdb_amount, pdb_diff_amount)
    rmsd_df = df.pivot(index=y_name, columns=x_name, values='rmsd')
    #print(rmsd_df)

    fig, ax = plt.subplots(figsize=(6,4))
    if xylabels:
        sns.heatmap(rmsd_df, cmap='viridis', cbar_kws={'label': r'RMSD ($\AA$)'}, 
                    xticklabels=xylabels, yticklabels=xylabels)
    else:
        sns.heatmap(rmsd_df, cmap='viridis', cbar_kws={'label': r'RMSD ($\AA$)'})
    #plt.imshow(rmsd_df, cmap='viridis')
    #plt.pcolor(rmsd_df, cmap='viridis')
    #plt.xlim(1,100)
    #plt.xticks(rmsd_df.columns)
    #ax.set_xticks([str(x) for x in rmsd_df.columns])
    #plt.ylim(0,100)
    #plt.yticks(np.arange(1, 5, 1))
    #plt.colorbar(label=r'RMSD ($\AA$)')
    
    # plt.tight_layout()
    # #plt.show()
    # plt.savefig('figures/we_dist_c2_clust_rmsd_1-75_39.png', dpi=300)

# compare km clusters to 2kod pdb paths (NMR ensemble)
# rmsd_list, pdb_comp_amount, pdb_diff = build_pdb_rmsd_matrix(nmr_pdb_paths, pdb_diff_path=km4_paths)
# plot_pdb_rmsd_matrix(rmsd_list, pdb_comp_amount, pdb_diff_amount=pdb_diff)

# compare monomer1 orientation to monomer2 orientations using 2KOD ensemble as reference point
#rmsd_list, pdb_comp_amount, pdb_diff = build_pdb_rmsd_matrix(nmr_pdb_paths, pdb_diff_path=mono2_paths)
#plot_pdb_rmsd_matrix(rmsd_list, pdb_comp_amount, pdb_diff_amount=pdb_diff)

# compare m2-m1 E175-W184 dist clusters of m01 2KOD 1us to 2KOD M01 and 1A43
#rmsd_list, pdb_comp_amount, pdb_diff = build_pdb_rmsd_matrix(main_pdb_paths, pdb_diff_path=m2_m1_m01_paths)
#plot_pdb_rmsd_matrix(rmsd_list, pdb_comp_amount, pdb_diff_amount=pdb_diff)

# compare C2 (2KOD style vector) clusters of m01 2KOD 1us to 2KOD M01 and 1A43
#rmsd_list, pdb_comp_amount, pdb_diff = build_pdb_rmsd_matrix(main_pdb_paths, pdb_diff_path=c2_m01_paths)
#plot_pdb_rmsd_matrix(rmsd_list, pdb_comp_amount, pdb_diff_amount=pdb_diff)

# compare  m2-m1 E175-W184 dist clusters and C2 (2KOD style vector) clusters of m01 2KOD 1us
#rmsd_list, pdb_comp_amount, pdb_diff = build_pdb_rmsd_matrix(m2_m1_m01_paths, pdb_diff_path=c2_m01_paths)
#plot_pdb_rmsd_matrix(rmsd_list, pdb_comp_amount, pdb_diff_amount=pdb_diff)

# complete comparison matrix
# total = main_pdb_paths + m2_m1_we_2kod_paths + m2_m1_we_1a43_paths+ c2_we_2kod_paths + c2_we_1a43_paths
# labels = ["2KOD", "2KOD Min", "1A43", "1A43 Min", 
#           "2KOD Dist C0", "2KOD Dist C1", "2KOD Dist C2", "1A43 Dist C0", "1A43 Dist C1", "1A43 Dist C2", 
#           "2KOD Angle C0", "2KOD Angle C1", "2KOD Angle C2", "1A43 Angle C0", "1A43 Angle C1", "1A43 Angle C2"
#           ]
# rmsd_list, pdb_comp_amount, pdb_diff = build_pdb_rmsd_matrix(total, pdb_diff_path=total)
# plot_pdb_rmsd_matrix(rmsd_list, pdb_comp_amount, pdb_diff_amount=pdb_diff, xylabels=labels)

plt.rcParams.update({'font.size': 14})
#plt.rcParams["figure.titleweight"] = "bold"
#plt.rcParams["font.weight"] = "bold"
#plt.rcParams["axes.labelweight"] = "bold"

# subset of comparison matrix
total = main_pdb_paths + c2_we_2kod_paths + c2_we_1a43_paths
labels = ["NMR: (2KOD)", "Xtal: (1A43)", 
          "NMR C-1", "NMR C-2", "NMR C-3", 
          "Xtal C-1", "Xtal C-2", "Xtal C-3"
          ]
rmsd_list, pdb_comp_amount, pdb_diff = build_pdb_rmsd_matrix(total, pdb_diff_path=total)
plot_pdb_rmsd_matrix(rmsd_list, pdb_comp_amount, pdb_diff_amount=pdb_diff, xylabels=labels)

plt.title("WE Data: Clustering of Helical Angle", fontweight="bold")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
#plt.show()
#plt.savefig('figures/we_c2_clust_rmsd_1-75_39.png', dpi=300, transparent=True)