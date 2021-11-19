"""
Extract the phi/psi or chi1/chi2 of a residue number of a pdb file.
Uses the Dihedral function of MDAnalysis with custom atom groups.
"""

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals

import pandas as pd

class Get_Dihedrals():
    """
    Take a pdb file or list of pdb files and calculate the dihedral angles.
    Should work for TRP, TYR, and PHE.
    """

    def __init__(self, pdb, resname):
        """
        Parameters
        ----------
        pdb : str
            Path to pdb file.
        resname : str
            Name of the residue to find dihedrals of; e.g. '4FW'.
        """
        # universe from mda of single pdb file
        self.u = mda.Universe(pdb, in_memory=True)
        # select all atoms with resname
        res = self.u.select_atoms(f"resname {resname}")
        # all unique 19F residue id(s) and segment id(s)
        self.res_ids = np.unique(res.resids)
        self.seg_ids = np.unique(res.segids)

    def calc_phi_psi(self):
        """
        Ramachandran.
        Phi = C(i-1), N, CA, C
        Psi = N, CA, C, N(i+1)

        Returns
        -------
        phipsi : ndarray
            2 column array where col1=phi and col2=psi values.
        phi : array
            Flat array of just the phi values.
        psi : array
            Flat array of just the psi values.
        """
        phi = []
        psi = []
        for seg in self.seg_ids:
            for res in self.res_ids:
                # when there are multiple atom locations from experiments, there can
                # be two possible locations for an atom, always pick the first location
                C_m1 = self.u.select_atoms(f"resid {res - 1} and name C and segid {seg}")[0]
                N = self.u.select_atoms(f"resid {res} and name N and segid {seg}")[0]
                CA = self.u.select_atoms(f"resid {res} and name CA and segid {seg}")[0]
                C = self.u.select_atoms(f"resid {res} and name C and segid {seg}")[0]
                N_p1 = self.u.select_atoms(f"resid {res + 1} and name N and segid {seg}")[0]

                phi.append(dihedrals.Dihedral([C_m1 + N + CA + C]).run().results.angles[0])
                psi.append(dihedrals.Dihedral([N + CA + C + N_p1]).run().results.angles[0])

        # stack by column so that col1=phi and col2=psi
        self.phipsi = np.hstack([phi, psi])

        # make into a flat 1D list of dihedral values
        self.phi = np.reshape(np.array(phi), -1)
        self.psi = np.reshape(np.array(psi), -1)

    def calc_chi1_chi2(self):
        """
        Janin.
        Chi1 = N, CA, CB, CG
        Chi2 = CA, CB, CG, CD1

        Returns
        -------
        chi1chi2 : ndarray
            2 column array where col1=chi1 and col2=chi2 values.
        chi1 : array
            Flat array of just the chi1 values.
        chi2 : array
            Flat array of just the chi2 values.
        """
        chi1 = []
        chi2 = []
        for seg in self.seg_ids:
            for res in self.res_ids:
                # when there are multiple atom locations from experiments, there can
                # be two possible locations for an atom, always pick the first location
                N = self.u.select_atoms(f"resid {res} and name N and segid {seg}")[0]
                CA = self.u.select_atoms(f"resid {res} and name CA and segid {seg}")[0]
                CB = self.u.select_atoms(f"resid {res} and name CB and segid {seg}")[0]
                CG = self.u.select_atoms(f"resid {res} and name CG and segid {seg}")[0]
                CD1 = self.u.select_atoms(f"resid {res} and name CD1 and segid {seg}")[0]

                chi1.append(dihedrals.Dihedral([N + CA + CB + CG]).run().results.angles[0])
                chi2.append(dihedrals.Dihedral([CA + CB + CG + CD1]).run().results.angles[0])

        # stack by column so that col1=chi1 and col2=chi2
        self.chi1chi2 = np.hstack([chi1, chi2])

        # make into a flat 1D list of dihedral values
        self.chi1 = np.reshape(np.array(chi1), -1)
        self.chi2 = np.reshape(np.array(chi2), -1)


def make_all_pdb_dh_dat_files():
    """
    Go through each 19F res class and generate then export both the stacked
    phi/psi and chi1/chi2 data arrays of all pdb files within the res class.
    """
    # key = 19F ipq res name, value = PDB identifier
    res_classes = {"w4f":"4FW", "w5f":"FTR", "w6f":"FT6", "w7f":"F7W",
                   "y3f":"YOF", "YDF":"F2Y", "f4f":"PFF", "ftf":"55I"}

    for key, val in res_classes.items():
        phipsi = np.array(["PHI", "PSI"])
        chi1chi2 = np.array(["CHI1", "CHI2"])
        parent = f"19F_pdbs/{key}"

        # calc for all pdb and ent files in directory
        for file in os.listdir(parent):
            if file.endswith((".pdb", ".ent")):
                dh = Get_Dihedrals(f"{parent}/{file}", val)
                dh.calc_phi_psi()
                dh.calc_chi1_chi2()
                phipsi = np.vstack((phipsi, dh.phipsi))
                chi1chi2 = np.vstack((chi1chi2, dh.chi1chi2))
    
        pd.DataFrame(phipsi).to_csv(f"19F_pdbs/dh_data/{key}_phipsi.tsv", sep="\t", header=False, index=False)
        pd.DataFrame(chi1chi2).to_csv(f"19F_pdbs/dh_data/{key}_chi1chi2.tsv", sep="\t", header=False, index=False)
        #np.savetxt(f"19F_pdbs/dh_data/{key}_phipsi.tsv", phipsi, delimiter="\t")
        #np.savetxt(f"19F_pdbs/dh_data/{key}_chi1chi2.tsv", chi1chi2, delimiter="\t")
                

# if __name__ == "__main__":
#     make_all_pdb_dh_dat_files() 