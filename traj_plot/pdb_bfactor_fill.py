"""
Python script to take input pdb file and write new with b-factor
replaced by another metric (rmsd, rmsf, csp, etc.) for 3-D visualization.
"""
import os

def replace_bf_pdb(pdb_input, pdb_output, data):
    """
    Replace column of b-factors (temperature factors) with data value.
    Data must be formatted: x = residue index, y = n digit float (TODO).
    """
    # read target coordinates
    pdb_coord = open(pdb_input, "r")
    pdb_coord_lines = pdb_coord.readlines()
    pdb_coord.close()

    temp = "temp_data_file.dat"
    # clean up data input - create new file (temp)
    data_input = open(temp, "w")
    for line in open(data).readlines():
        # skip lines that begin with "#" or is blank or blank with only tab
        # TODO: also skip lines that begin with "#" as the first non-empty char in the line
        if line[:1] == "#" or line == "\n" or line == "\t\n":
            continue
        else:
            data_input.write(line)
    data_input.close()

    # list comprehension to import columns of tsv
    x = [float(((x.split())[0])) for x in open(temp).readlines()]
    y = [str((y.split())[1])[0:5] for y in open(temp).readlines()] # read 5 char

    # clean up temp file
    os.remove(temp)

    # make dictionary of res_num (1-176) (x) and y value
    data_zip = dict(zip(x, y))
    #print(data_zip)

    # create new coordinate file
    new_pdb = open(pdb_output, "w")

    # go through each line in target pdb file
    for line in pdb_coord_lines:

        # focus on lines that have atomic information
        if line.startswith('ATOM'):

            chain = str(line[21])
            res_num = float(line[22:26])
            #occupancy = float([54:60])
            b_fact = str(line[60:66])

            # TODO: update this for Xtal pdb files, replaces original bf value with 0.00
            # if b_fact != "  0.00":
            #     new_line = line.replace(b_fact, "  0.00")

            # adjust residue naming conventions to be 1-176 (same as data value indicies)
            if chain == "A":
                res_num = res_num - 143
            if chain == "B":
                res_num = res_num - 143 + 88

            # replace b factor with data for matching res num
            for res, value in data_zip.items():
                if res_num == res:
                    # data "value" should be 5 char, with leading " " = 6 char = b-factor field len
                    # if the data value is less than 5 char, add trailing spaces
                    new_line = line.replace(b_fact, " " + str(value) + ((5 - len(value)) * " "))

            # if the current residue of pdb file has a cooresponding residue value in data file
            # this ensures that datasets with missing data for any residues will be skipped (bf remains 0.00)
            if res_num in data_zip.keys():
                new_pdb.write(new_line)

            # also write out 'ATOM' lines without matching data present
            else:
                new_pdb.write(line)
            
        # write all other non 'ATOM' lines as well
        else:
            new_pdb.write(line)

pdb_file = "/Users/Darian/Research/Projects/hiv1_capsid_sim/pdb_files/2KOD/m01/m01_2kod_pH_amb.pdb"
#replace_bf_pdb(pdb_file, "test_out.pdb", "m01/m01_res_rmsf_copy.dat")