import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align, rms
import numpy as np
import matplotlib.pyplot as plt

# choose your 2kod model number
m = 'm04'

# load parameters and trajectories
m01 = mda.Universe(f'2kod_{m}/{m}_2kod_dry.prmtop', f'2kod_{m}/{m}_2kod_1us_10i_dry.ncdf', \
    in_memory=True, in_memory_step=100, verbose=True)


def single_traj_pw_rmsd(m):
    """
    Single traj test of MDAnalysis code.
    """
    # load parameters and trajectories
    traj = mda.Universe(f'2kod_m0{m}/m0{m}_2kod_dry.prmtop', f'2kod_m0{m}/m0{m}_2kod_1us_10i_dry.ncdf', \
        in_memory=True, in_memory_step=100, verbose=True)
    
    # only look at non-termini
    #traj = traj.select_atoms('resnum 6:75 and not name H*', 'resnum 94:163 and not name H*')

    # align traj to minimize rmsd, should be alligned from cpptraj
    align.AlignTraj(traj, traj, select=('resnum 6:75 and not name H*', 'resnum 94:163 and not name H*'), in_memory=True).run()

    # calc pairwise rmsd matrix
    matrix = diffusionmap.DistanceMatrix(traj, select='resnum 6:75 or resnum 94:163 and not name H*').run()
    #print(matrix.dist_matrix.shape)

    # plot a heatmap
    plt.figure(m)
    plt.imshow(matrix.dist_matrix, vmin=0, vmax=10, cmap='viridis')
    plt.xlabel('Time (ns)')
    #plt.xlim(0,100)
    #plt.xticks(np.arange(0, 101, 1))
    #plt.ylim(0,100)
    #plt.yticks(np.arange(0, 101, 1))
    plt.ylabel('Time (ns)')
    plt.colorbar(label=r'RMSD ($\AA$)')

    plt.savefig(f'm0{m}_1us_pw_rmsd.png', dpi=300)


for i in range(1,5):
    single_traj_pw_rmsd(i)



def two_traj_pw_rmsd(traj_1, traj_2):
    """ 
    Multi traj MDAnalysis test code.
    """
    # empty 2-D array to become trajectory comparison results
    prmsd = np.zeros((len(traj_1.trajectory),  # y-axis
                    len(traj_2.trajectory)))   # x-axis

    # iterate each frame of traj_1 and compare to traj_2 and store
    for i, frame_open in enumerate(traj_1.trajectory):
        r = rms.RMSD(traj_2, traj_1, select='name CA',
                    ref_frame=i).run()
        prmsd[i] = r.rmsd[:, -1]  # select 3rd column with RMSD values

    # plot heatmap    
    plt.imshow(prmsd, cmap='viridis')
    plt.xlabel('Frame (adk_closed)')
    plt.ylabel('Frame (adk_open)')
    plt.colorbar(label=r'RMSD ($\AA$)')

    plt.savefig("figures/two_traj_test.png")
    
#two_traj_pw_rmsd(m01, m02)
