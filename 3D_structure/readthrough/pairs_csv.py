from argparse import ArgumentParser, RawTextHelpFormatter
import pandas as pd
import os
import numpy as np
from Bio import SeqIO, pairwise2
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBParser, Superimposer
import warnings
import json
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
warnings.filterwarnings('ignore')



############### functions ###############

#### RMSD ####

    
# to calculate rmsd using Bio.PDB
# predicted_pdb: moving pdb comparing to reference, need to have the same number of CA atoms as the ref_pdb
def rmsd_pdb(predicted_pdb, ref_pdb):
    parser = PDBParser()
    struct_ref = parser.get_structure(os.path.basename(ref_pdb), ref_pdb)
    struct_predicted = parser.get_structure(os.path.basename(predicted_pdb), predicted_pdb)
    fixed = [atom for atom in struct_ref[0].get_atoms() if atom.name == "CA"]
    moving = [atom for atom in struct_predicted[0].get_atoms() if atom.name == "CA"]
    sup = Superimposer()
    # sets the fixed and moving atom lists
    # finds the rotation and translation matrix that best superimposes the moving atoms onto fixed atoms
    sup.set_atoms(fixed, moving)
    # applies the calculated rotation and translation to all atoms in the second structure
    # superimposing it onto the first structure
    #sup.apply(struct_predicted[0].get_atoms())
    sup.apply(moving)

    return sup.rms


def rmsd_point(coordinate1, coordinate2):
    # 3D coordinates, example: array([ 11.925492,  10.070204, -12.518902], dtype=float32)
    # this is a function to calculate rmsd for a single point
    x1 = coordinate1[0]
    y1 = coordinate1[1]
    z1 = coordinate1[2]
    x2 = coordinate2[0]
    y2 = coordinate2[1]
    z2 = coordinate2[2]
    value = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    return value


def rmsd_list(coordinates1, coordinates2):
    # list of 3D coordinates, should have the same number of coordinates
    # this is a function to calculate rmsd for two list of 3D coordinates
    length = len(coordinates1)
    values = []
    for i in range(length):
        x1 = coordinates1[i][0]
        y1 = coordinates1[i][1]
        z1 = coordinates1[i][2]
        x2 = coordinates2[i][0]
        y2 = coordinates2[i][1]
        z2 = coordinates2[i][2]
        value = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        values.append(value)
    rmsd = np.sqrt(sum(values)/length)
    return rmsd



#### Pairwise Sequence Alignment ####
    

# load fasta file to get the sequence
def load_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        record = next(SeqIO.parse(f, 'fasta'))
        return record.seq


# extract common residues from two sequences
def common_residues(seq1, seq2):
    # inputs are aligned sequences, including gaps, length should be the same
    # outputs are the common residues of the two sequences, indices are the original ones (without gaps)
    common_indices1 = []
    common_indices2 = []
    res_idx1 = 0
    res_idx2 = 0
    for a, b in zip(seq1, seq2):
        if a == b:
            res_idx1 += 1
            res_idx2 += 1
            common_indices1.append(res_idx1)
            common_indices2.append(res_idx2)
        elif (a != b) and (a != '-') and (b != '-'):
            res_idx1 += 1
            res_idx2 += 1
        elif (a == '-') and (b != '-'):
            res_idx2 += 1
        elif (a != '-') and (b == '-'):
            res_idx1 += 1
        else:
            print(a, b)
    if len(common_indices1) != len(common_indices2):
        print("Two indices have different length!")
        return 1
    return common_indices1, common_indices2

    
# return common indices for each sequence
def pairwise_SA(moving_fasta, ref_fasta):
    # load sequences from fasta files
    seq_moving = load_fasta(moving_fasta)
    seq_ref = load_fasta(ref_fasta)

    # pairwise sequence alignment
    alignments = pairwise2.align.globalxx(seq_moving, seq_ref)
    aligned_seq_moving, aligned_seq_ref = alignments[0][:2]

    # extract common residues
    indices_moving, indices_ref = common_residues(aligned_seq_moving, aligned_seq_ref)

    return indices_moving, indices_ref


#### read residue types from pairs ####
def residue_type(moving_fasta, ref_fasta):
    # load sequences from fasta files
    seq_moving = load_fasta(moving_fasta)
    seq_ref = load_fasta(ref_fasta)

    # pairwise sequence alignment
    alignments = pairwise2.align.globalxx(seq_moving, seq_ref)
    aligned_seq_moving, aligned_seq_ref = alignments[0][:2]

    # extract common residues
    indices_moving, indices_ref = common_residues(aligned_seq_moving, aligned_seq_ref)

    # initialize dictionary to save residue types
    residue_moving = {}
    residue_ref = {}

    # save the residue
    for i in range(len(indices_ref)):
        idx_moving = indices_moving[i] 
        idx_ref = indices_ref[i]
        residue_moving[idx_ref] = seq_moving[idx_moving-1] # use reference index
        residue_ref[idx_ref] = seq_ref[idx_ref-1] # change 1-based index to 0-based index

    # check if the two dictionaries are the same
    if residue_moving != residue_ref:
        print("residue inconsistency!")
        return 1

    # only need to return one dictionary
    return residue_ref
    
    
#### RMSD with pair ####
    
# to calculate rmsd per residue using Bio.PDB
def rmsd_pdb_perResidue(moving_pdb, ref_pdb, moving_fasta, ref_fasta):
    # load PDB structures
    parser = PDBParser()
    struct_moving = parser.get_structure(os.path.basename(moving_pdb), moving_pdb)
    struct_ref = parser.get_structure(os.path.basename(ref_pdb), ref_pdb)

    # extract common residues
    indices_moving, indices_ref = pairwise_SA(moving_fasta, ref_fasta)
    
    # get CA atoms from PDBs
    moving = [atom for atom in struct_moving[0].get_atoms() if (atom.name == "CA") and (atom.full_id[3][1] in indices_moving)]
    fixed = [atom for atom in struct_ref[0].get_atoms() if (atom.name == "CA") and (atom.full_id[3][1] in indices_ref)]

    # check if the two structures have the same number of length
    if len(moving) != len(fixed):
        print("Two structures have different numbers of residues!")
    
    # get the fixed coordinates
    coords_fixed = []
    for i in range(len(fixed)):
        coords_fixed.append(fixed[i].get_coord())
    # get the moving coordinates
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    sup.apply(moving)
    coords_moving = []
    for i in range(len(moving)):
        coords_moving.append(moving[i].get_coord())

    # calculate rmsd per residue (CA)
    rmsd_perResidue = {}
    for i in range(len(coords_fixed)):
        residue_id = fixed[i].full_id[3][1]
        rmsd_perResidue[residue_id] = rmsd_point(coords_fixed[i], coords_moving[i])

    return rmsd_perResidue


# to calculate rmsd overall using Bio.PDB
# return rmsd between two structures
# indices consistent with the previous function "rmsd_pdb_perResidue"
def rmsd_pdb_overall(moving_pdb, ref_pdb, moving_fasta, ref_fasta):
    # load PDB structures
    parser = PDBParser()
    struct_moving = parser.get_structure(os.path.basename(moving_pdb), moving_pdb)
    struct_ref = parser.get_structure(os.path.basename(ref_pdb), ref_pdb)

    # extract common residues
    indices_moving, indices_ref = pairwise_SA(moving_fasta, ref_fasta)
    
    # get CA atoms from PDBs
    moving = [atom for atom in struct_moving[0].get_atoms() if (atom.name == "CA") and (atom.full_id[3][1] in indices_moving)]
    fixed = [atom for atom in struct_ref[0].get_atoms() if (atom.name == "CA") and (atom.full_id[3][1] in indices_ref)]

    # check if the two structures have the same number of length
    if len(moving) != len(fixed):
        print("Two structures have different numbers of residues!")
    
    # get the fixed coordinates
    coords_fixed = []
    for i in range(len(fixed)):
        coords_fixed.append(fixed[i].get_coord())
    # get the moving coordinates
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    sup.apply(moving)
    coords_moving = []
    for i in range(len(moving)):
        coords_moving.append(moving[i].get_coord())

    rmsd = rmsd_list(coords_moving, coords_fixed)

    return rmsd



#### delta delta G ####
    
# extract delta G per residue out of sc file
def dG_perResidue(sc_path):
    dG_perResidue = {}
    with open(sc_path, 'r') as f:
        for count, line in enumerate(f.readlines()):
            if(count != 0):
                line = line.strip("\n").split()
                id = int(line[23].split("_")[1])
                score = float(line[22])
                dG_perResidue[id] = score
    return dG_perResidue

# to calculate ddG within pairs
def ddG_perResidue(moving_sc, ref_sc, moving_fasta, ref_fasta):
    #### input csv files !!! ####

    # extract common residues
    indices_moving, indices_ref = pairwise_SA(moving_fasta, ref_fasta)

    # read score file
    dG_moving = pd.read_csv(moving_sc)
    dG_ref = pd.read_csv(ref_sc)

    # initiate dG and ddG dictionary
    dG_iso_perResidue = {}
    dG_ref_perResidue = {}
    ddG_perResidue = {}
    
    # dG and ddG
    for i in range(len(indices_ref)):
        idx_moving = indices_moving[i]
        score_moving = dG_moving[dG_moving['residue_id']==idx_moving]['dG'].iloc[0]
        idx_ref = indices_ref[i]
        score_ref = dG_ref[dG_ref['residue_id']==idx_ref]['dG'].iloc[0]
        ddG = round(score_moving - score_ref, 3)

        # use index for reference as residue id
        dG_iso_perResidue[idx_ref] = score_moving
        dG_ref_perResidue[idx_ref] = score_ref
        ddG_perResidue[idx_ref] = ddG

    return ddG_perResidue, dG_iso_perResidue, dG_ref_perResidue
    

#### extract fa_sol term out of sc files ####
def fa_sol_perResidue(sc_path):
    fa_sol_perResidue = {}
    with open(sc_path, 'r') as f:
        for count, line in enumerate(f.readlines()):
            if (count != 0):
                line = line.strip("\n").split()
                id = int(line[23].split("_")[1])
                fa_sol = float(line[5])
                fa_sol_perResidue[id] = fa_sol
    return fa_sol_perResidue

    
# get fa_sol within pairs
def get_pair_fa_sol(moving_sc, ref_sc, moving_fasta, ref_fasta):
    #### input csv files !!! ####

    # extract common residues
    indices_moving, indices_ref = pairwise_SA(moving_fasta, ref_fasta)

    # read score file
    fa_sol_moving = pd.read_csv(moving_sc)
    fa_sol_ref = pd.read_csv(ref_sc)

    # initiate fa_sol dictionary
    fa_sol_iso_perResidue = {}
    fa_sol_ref_perResidue = {}

    # fa_sol
    for i in range(len(indices_ref)):
        idx_moving = indices_moving[i]
        score_moving = fa_sol_moving[fa_sol_moving['residue_id']==idx_moving]['fa_sol'].iloc[0]
        idx_ref = indices_ref[i]
        score_ref = fa_sol_ref[fa_sol_ref['residue_id']==idx_ref]['fa_sol'].iloc[0]

        # use index for reference as residue id
        fa_sol_iso_perResidue[idx_ref] = score_moving
        fa_sol_ref_perResidue[idx_ref] = score_ref

    return fa_sol_iso_perResidue, fa_sol_ref_perResidue

    
# get tags with the lowest relax score
# here it generates a tag file (basically a txt file)
# that contains 3 ids with lowest score by default

def get_lowestTag(sc_file, tag_file, num = 3):
    if sc_file==tag_file:
        print("score file and tag file should be different!")
        return 1
    # sc file is the one generated by Rosetta Relax
    # tag file is a path, not a folder, need to specify the file name
    scores_and_ids = pd.DataFrame(columns = ['score', 'id'])
    with open(sc_file, "r") as f:
        for count, line in enumerate(f.readlines()):
            if (count != 0) and (count != 1): # excluded the first two lines
                line = line.strip("\n")
                line = line.split()
                scores_and_ids.loc[len(scores_and_ids)] = [float(line[1]), line[23]]
    scores_and_ids = scores_and_ids.sort_values(by = 'score', ascending = True)
    lowest_ids = scores_and_ids['id'].head(num)
    with open(tag_file, 'w') as f:
        for id in lowest_ids:
            f.write(f"{id}\n")
    return 0


# get the lowest relax score
def get_lowestScore(sc_file):
    # very similar to the function get_lowestTag
    # sc file is the one generated by Rosetta Relax
    scores_and_ids = pd.DataFrame(columns = ['score', 'id'])
    with open(sc_file, 'r') as f:
        for count, line in enumerate(f.readlines()):
            if (count != 0) and (count != 1): # excluded the first two lines
                line = line.strip("\n")
                line = line.split()
                scores_and_ids.loc[len(scores_and_ids)] = [float(line[1]), line[23]]
    scores_and_ids = scores_and_ids.sort_values(by = 'score', ascending = True)
    lowest_score = float(scores_and_ids['score'].head(1))
    lowest_tag = scores_and_ids['id'].head(1).iloc[0]
    return lowest_score, lowest_tag
    

#### plddt score (with pair) ####
    
# get plddt score from "ranking_debug.json"
def get_plddt(af_ranking_file):
    af_ranking = json.load(open(af_ranking_file))
    ave_plddt = format(sum(af_ranking['plddts'].values()) / len(af_ranking['plddts'].values()), '.3f')
    return float(ave_plddt)


# get lowest plddt score from "ranking_debug.json"
def get_plddt_highest(af_ranking_file):
    af_ranking = json.load(open(af_ranking_file))
    plddt = max(af_ranking['plddts'].values())
    return plddt

# get plddt per residue
# based on the original residue id within the corresponding pdb file
def get_plddt_perResidue(result_model_pkl):
    plddt_dic = {} # create a dictionary to store the plddt score
    plddt = pd.read_pickle(result_model_pkl)
    plddt = plddt['plddt']
    for i in range(len(plddt)):
        plddt_dic[i+1] = plddt[i]
    return plddt_dic


# get model name for highest ranking
def get_ranked_0_model(af_ranking_file):
    af_ranking = json.load(open(af_ranking_file))
    model_name = af_ranking['order'][0]
    return model_name


# get plddt per residue, index based on pairwise sequence alignment
def plddt_pair_perResidue(moving_pkl, ref_pkl, moving_fasta, ref_fasta):
    ### input pkl files

    # extract common residues
    indices_moving, indices_ref = pairwise_SA(moving_fasta, ref_fasta)

    # read pkl files
    moving_plddt = get_plddt_perResidue(moving_pkl)
    ref_plddt = get_plddt_perResidue(ref_pkl)

    # initiate plddt dictionary
    moving_plddt_perResidue = {}
    ref_plddt_perResidue = {}
    
    # change the index based on sequence alignment
    for i in range(len(indices_ref)):
        idx_moving = indices_moving[i]
        idx_ref = indices_ref[i]
        plddt_score_moving = moving_plddt[idx_moving]
        plddt_score_ref = ref_plddt[idx_ref] # get plddt score from dictionary
        moving_plddt_perResidue[idx_ref] = plddt_score_moving # use reference index for consistency
        ref_plddt_perResidue[idx_ref] = plddt_score_ref

    return moving_plddt_perResidue, ref_plddt_perResidue




    
#### RSA from DSSP ####


# get relative accessible surface
def get_dssp(pdb_file):
    # need to import #
    # from Bio.PDB.DSSP import DSSP
    # from Bio.PDB import PDBParser
    p = PDBParser()
    structure = p.get_structure(str(pdb_file), pdb_file)
    model = structure[0]

    # count residues
    #n = 0
    #for res in model.get_residues():
        #n += 1
    
    dssp = DSSP(model, pdb_file, dssp = 'mkdssp')
    #if len(dssp.keys()) != n:
        #print(os.path.basename(pdb_file) + " dssp length different from pdb!")
        #return 1
    return dssp


# change from residue index to dssp index
def fix_index(pdb_file, chain_id, residue_id):
    # from Bio.PDB.DSSP import dssp_dict_from_pdb_file

    # define chain
    #chain_id = os.path.basename(pdb_file).split('.')[0].split('_')[1]

    # get DSSP index using residue index
    dssp_tuple = dssp_dict_from_pdb_file(pdb_file)
    dssp_dict = dssp_tuple[0]
    key = (chain_id, (' ', residue_id, ' '))
    if key in dssp_dict.keys():
        dssp_id = dssp_dict[chain_id, (' ', residue_id, ' ')][5]
        # reference: https://github.com/biopython/biopython/blob/master/Bio/PDB/DSSP.py
        return dssp_id
    return None

    
# get single RSA value given the residue ID
def get_single_rsa(pdb_file, chain_id, residue_id):
    dssp_data = get_dssp(pdb_file)
    dssp_id = fix_index(pdb_file, chain_id, residue_id)

    if dssp_id == None:
        return None

    for key in dssp_data.keys():
        if dssp_data[key][0] == dssp_id:
            # dssp_data[key][0] is DSSP index
            # according to https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html
            return dssp_data[key][3]
    return None

    
# get RSA per residue, index based on pairwise sequence alignment
def get_pair_rsa(moving_pdb, ref_pdb, moving_fasta, ref_fasta, moving_chain_id, ref_chain_id):

    # need to import #
    # from Bio import SeqIO, pairwise2
    # from Bio.PDB.DSSP import DSSP
    # from Bio.PDB import PDBParser

    # extract common residues
    indices_moving, indices_ref = pairwise_SA(moving_fasta, ref_fasta)

    # get dssp
    #moving_dssp = get_dssp(moving_pdb)
    #ref_dssp = get_dssp(ref_pdb)

    # initiate rsa dictionary
    moving_rsa_perResidue = {}
    ref_rsa_perResidue = {}

    # index based on reference sequence
    for i in range(len(indices_ref)):
        idx_moving = indices_moving[i]
        idx_ref = indices_ref[i]
        
        # get RSA
        rsa_moving = get_single_rsa(moving_pdb, moving_chain_id, idx_moving)
        if rsa_moving != None:
            moving_rsa_perResidue[idx_ref] = rsa_moving # save in dictionary, use reference index for consistency
            
        rsa_ref = get_single_rsa(ref_pdb, ref_chain_id, idx_ref)
        if rsa_ref != None:
            ref_rsa_perResidue[idx_ref] = rsa_ref

    return moving_rsa_perResidue, ref_rsa_perResidue


############################################



parser = ArgumentParser(description = "Stability Metrics for Each Pair, arguments",
                        formatter_class = RawTextHelpFormatter)

# Support level input
parser.add_argument('--sl', type = str, required = True,
                    help = "Support Level (Proteomics_Evidence, RPF_Evidence, No_Evidence)")

# Protein ID input
parser.add_argument('--pid', type = str, required = True,
                    help = "Protein ID, [ProteinName]_[Index]")

options = parser.parse_args()


# Path storing all AlphaFold and RosettaRelax results
res_path = "/scicore/home/zavolan/zhu0006/3D_structure/structural_stability_of_novel_orfs_hcc/structure_res"
# Current support level result path
res_path = os.path.join(res_path, options.sl)

# Path storing all fasta files
fasta_path = "/scicore/home/zavolan/zhu0006/3D_structure/structural_stability_of_novel_orfs_hcc/fasta"
# Current support level fasta path
fasta_path = os.path.join(fasta_path, options.sl)

# Some file names
tag_filename = "lowest.tag"
sc_csv_filename = "relax_lowestScore_perRes.csv"

# Path to save csv results
pairs_csv_path = "/scicore/home/zavolan/zhu0006/3D_structure/structural_stability_of_novel_orfs_hcc/pairs_csv"


####### Start calculating #######

# get ids
id_C = options.pid + "_C"
id_R = options.pid + "_R"

# alphafold results path
res_path_C = os.path.join(res_path, id_C)
res_path_R = os.path.join(res_path, id_R)

# fasta paths
fasta_name_C = id_C + ".fasta"
fasta_path_C = os.path.join(fasta_path, fasta_name_C)
fasta_name_R = id_R + ".fasta"
fasta_path_R = os.path.join(fasta_path, fasta_name_R)

# tag paths (with lowest relax score)
tag_path_C = os.path.join(res_path_C, tag_filename)
tag_path_R = os.path.join(res_path_R, tag_filename)

# structure basenames (with lowest relax score)
with open(tag_path_C, 'r') as f:
    relax_basename_C = f.readline().strip()
with open(tag_path_R, 'r') as f:
    relax_basename_R = f.readline().strip()

# pdb file paths (with lowest relax score)
pdb_name_C = relax_basename_C + ".pdb"
pdb_path_C = os.path.join(res_path_C, pdb_name_C)
pdb_name_R = relax_basename_R + ".pdb"
pdb_path_R = os.path.join(res_path_R, pdb_name_R)

# score per residue csv file paths
sc_path_C = os.path.join(res_path_C, sc_csv_filename)
sc_path_R = os.path.join(res_path_R, sc_csv_filename)

#### indices for isoform ####
indices_R, indices_C = pairwise_SA(moving_fasta=fasta_path_R, ref_fasta=fasta_path_C) # return moving indices first (readthrough), then reference indices (canonical)
indices_dict = dict(zip(indices_C, indices_R)) # use canonical indices as keys
# initialize the dataframe
pair_df = pd.DataFrame({
    'Residue_C': indices_dict.keys(),
    'Residue_R': indices_dict.values()
})

#### add residue types ####
residues = residue_type(moving_fasta=fasta_path_R, ref_fasta=fasta_path_C)
# combine into the dataframe
pair_df['ResidueType'] = pair_df['Residue_C'].map(residues)

#### RMSD ####
rmsd = rmsd_pdb_perResidue(moving_pdb=pdb_path_R, ref_pdb=pdb_path_C, moving_fasta=fasta_path_R, ref_fasta=fasta_path_C)
# combine into the dataframe
pair_df['RMSD'] = pair_df['Residue_C'].map(rmsd)

#### dG and ddG ####
ddG, dG_R, dG_C = ddG_perResidue(moving_sc=sc_path_R, ref_sc=sc_path_C, moving_fasta=fasta_path_R, ref_fasta=fasta_path_C) #ddG=moving-ref
# combine into the dataframe
pair_df['ddG'] = pair_df['Residue_C'].map(ddG)
pair_df['dG_R'] = pair_df['Residue_C'].map(dG_R)
pair_df['dG_C'] = pair_df['Residue_C'].map(dG_C)

#### plddt per residue ####
model_name_C = get_ranked_0_model(af_ranking_file=os.path.join(res_path_C, "ranking_debug.json")) # get original model name with the highest confidence score
pkl_name_C = "result_" + model_name_C + ".pkl"
pkl_path_C = os.path.join(res_path_C, pkl_name_C)
model_name_R = get_ranked_0_model(af_ranking_file=os.path.join(res_path_R, "ranking_debug.json"))
pkl_name_R = "result_" + model_name_R + ".pkl"
pkl_path_R = os.path.join(res_path_R, pkl_name_R)
plddt_R, plddt_C = plddt_pair_perResidue(moving_pkl=pkl_path_R, ref_pkl=pkl_path_C, moving_fasta=fasta_path_R, ref_fasta=fasta_path_C)
# combine into the dataframe
pair_df['plddt_R'] = pair_df['Residue_C'].map(plddt_R)
pair_df['plddt_C'] = pair_df['Residue_C'].map(plddt_C)

#### RSA per residue ####
rsa_R, rsa_C = get_pair_rsa(moving_pdb=pdb_path_R, ref_pdb=pdb_path_C,
                            moving_fasta=fasta_path_R, ref_fasta=fasta_path_C,
                            moving_chain_id="A", ref_chain_id="A")
# combine into the dataframe
pair_df['rsa_R'] = pair_df['Residue_C'].map(rsa_R)
pair_df['rsa_C'] = pair_df['Residue_C'].map(rsa_C)

# print process
print(f"Processed {id_C} and {id_R}")

# save to csv
csv_name = options.pid + ".csv"
csv_path = os.path.join(pairs_csv_path, options.sl, csv_name)
pair_df.to_csv(csv_path, index = False)