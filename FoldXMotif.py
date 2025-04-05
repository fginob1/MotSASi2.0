#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from Bio.PDB import PDBParser
import re
import argparse
import gzip
from Bio.SeqUtils import seq1
from tqdm import tqdm
import statistics
import statsmodels.formula.api as smf
import statsmodels.api as sm

class pdb_match_object(object):
    """
    We think the pair motif-receptor domain as an object with attributes reflecting this information
    The pdb_match object carry all this valuable information that we are going to use in the future
    when running FoldX
    """
    def __init__(self, pdb, chain, peptide_start, hit, structure, info_gaps, pdb_chain):
        self.PDB = pdb
        self.chain = chain
        # we use this strategy to account for gaps for preventing missassigning of motif position numbering
        dif_total = 0
        for info_gap in info_gaps:
            print(info_gap)
            #if hit.start() < info_gap["Start"]:
                #continue
            #elif peptide_start+hit.start() > info_gap["End"]:
                #dif_total += info_gap["Difference"]-1
            if peptide_start+hit.start() > info_gap["Start"]:
                dif_total += info_gap["Difference"]-1
        print(peptide_start, hit.start(), dif_total)
        self.start_position = peptide_start + hit.start() + dif_total
        self.motif_seq = hit.group()
        self.motif_pdb = pdb_chain[hit.start():hit.end()]
        # this is a strategy to determine the chain of the receptor domain, basically we determine the nearest
        # chain to the residue in the middle of the motif (as alpha carbon)
        min_distance = 1000000000000000000000000
        middle_point = self.start_position + (len(self.motif_seq)//2)
        for model1 in structure:
            for chain1 in model1:
                if chain1.id == chain:
                    for residue1 in chain1:
                        if residue1.id[1] == middle_point:
                            for model2 in structure:
                                for chain2 in model2:
                                    if chain2.id != self.chain:
                                        for residue2 in chain2:
                                            if residue2.id[0] == " " and len(residue2.get_resname().strip()) == 3:
                                                for atom2 in residue2:
                                                    if atom2.fullname.strip() == "CA":
                                                        distance = residue1["CA"]-residue2["CA"]
                                                        if distance < min_distance:
                                                            min_distance = distance
                                                            self.partner_chain = chain2.id

    def __eq__(self, other) :
        return self.__dict__ == other.__dict__


def pdb_foldx_parameters(match, motif_re, motif_name, pdb_db_path):
    """
    Given a true positive motif object, get the pdb_match objects that define the interaction to drive
    a binding free energy analysis
    In: true positive motif object, motif name and re, pdb database path
    Out: pdb_match objects list
    """

    # list where storing pdb_match objects
    pdb_matchs = []

    # iterate through crystal structures, parsing them
    for pdb in match.PDBs:
        parser = PDBParser(PERMISSIVE=1)
        structure_id = pdb
        filename = f"{pdb_db_path}pdb{pdb.lower()}.ent.gz"
        structure = parser.get_structure(structure_id, gzip.open(filename, "rt"))
        chains = []
        for model in structure:
            # when iterating through chains, we are going to store residues on them in two ways:
            # firstly, as defined by the standard 20 aminoacids
            # secondly, incorporating alternative residues that can appear in a crystal, e.g. phosphorylated residues
            # if gaps are detected, we are going to retain that info
            for chain in model:
                residues_seq = []
                residues_pdb = []
                info_gaps = []
                for residue in chain:
                    # checking for standard and special residues while iterating and keeping them in the lists
                    # if it is the first residue, we keep the residue crystal structure numberation and skip gap testing
                    # always keep the last crystal structure residue parsed number
                    if residue.id[0] in [" ", "H_PTR", "H_TPO", "H_SEP", "H_MSE"] and len(residue.get_resname().strip()) == 3:
                        # phosphorylated tyrosine
                        if residue.get_resname() == "PTR":
                            residues_pdb.append("y")
                            residues_seq.append("Y")
                        # phosphorylated threonine
                        elif residue.get_resname() == "TPO":
                            residues_pdb.append("p")
                            residues_seq.append("T")
                        # phosphorylated serine
                        elif residue.get_resname() == "SEP":
                            residues_pdb.append("s")
                            residues_seq.append("S")
                        # seleno-methionine
                        elif residue.get_resname() == "MSE":
                            # no matter foldx, keep a methionine
                            residues_pdb.append("M")
                            residues_seq.append("M")
                        else:
                            residues_pdb.append(seq1(residue.get_resname()))
                            residues_seq.append(seq1(residue.get_resname()))
                        if len(residues_seq) == 1:
                            peptide_start = residue.id[1]
                            previous_aa = residue.id[1]
                            continue

                    # info_gap is a dictionary with the gap coordinates
                    if residue.id[1] != previous_aa + 1 and residue.id[0] in [" ", "H_PTR", "H_TPO", "H_SEP", "H_MSE"] and len(residue.get_resname().strip()) == 3:
                        info_gap = {}
                        dif = residue.id[1] - previous_aa
                        info_gap["Start"] = previous_aa
                        info_gap["End"] = residue.id[1]
                        info_gap["Difference"] = dif
                        info_gaps.append(info_gap)

                    previous_aa = residue.id[1]

                # as motifs are generally present in small chains, we look for them
                if len(residues_seq) < 50:
                    seq_chain = "".join(residues_seq)
                    pdb_chain = "".join(residues_pdb)
                    # look for matches in the sequence chain, where residues are as the standard 20,
                    # and then we build a pdb_match object with it and keep it
                    for hit in re.finditer(motif_re, str(seq_chain)):
                        pdb_match = pdb_match_object(pdb, chain.id, peptide_start, hit, structure, info_gaps, pdb_chain)
                        if pdb_match not in pdb_matchs:
                            pdb_matchs.append(pdb_match)
                elif motif_name == "LIG_ARL_BART_1" and "3DOE" in match.PDBs and chain.id == "A":
                    seq_chain = "".join(residues_seq)
                    pdb_chain = "".join(residues_pdb)
                    # look for matches in the sequence chain, where residues are as the standard 20,
                    # and then we build a pdb_match object with it and keep it
                    for hit in re.finditer(motif_re, str(seq_chain)):
                        pdb_match = pdb_match_object(pdb, chain.id, peptide_start, hit, structure, info_gaps, pdb_chain)
                        if pdb_match not in pdb_matchs:
                            pdb_matchs.append(pdb_match)

    return pdb_matchs

def alphafold_foldx_parameters(motif_re, motif_name, row, cwd_path):
    """
    Given a true positive motif object, get the pdb_match objects that define the interaction to drive
    a binding free energy analysis
    In: true positive motif object, motif name and re, pdb database path
    Out: pdb_match objects list
    """

    # iterate through crystal structures, parsing them
    parser = PDBParser(PERMISSIVE=1)
    structure_id = "alphafold"
    ranking = row["rank"]
    for file in os.listdir(f"{cwd_path}Motifs/{motif_name}/Seed/"):
        if file.startswith(motif_name) and file.endswith(f"AF2_rank_00{ranking}_model_{row.model}.pdb"):
            filename = f"{cwd_path}Motifs/{motif_name}/Seed/{file}"
            break

    structure = parser.get_structure(structure_id, filename)
    chains = []
    for model in structure:
        # when iterating through chains, we are going to store residues on them in two ways:
        # firstly, as defined by the standard 20 aminoacids
        # secondly, incorporating alternative residues that can appear in a crystal, e.g. phosphorylated residues
        # if gaps are detected, we are going to retain that info
        for chain in model:
            residues_seq = []
            residues_pdb = []
            info_gaps = []
            for residue in chain:
                # checking for standard and special residues while iterating and keeping them in the lists
                # if it is the first residue, we keep the residue crystal structure numberation and skip gap testing
                # always keep the last crystal structure residue parsed number
                if residue.id[0] in [" ", "H_PTR", "H_TPO", "H_SEP", "H_MSE"] and len(residue.get_resname().strip()) == 3:
                    # phosphorylated tyrosine
                    if residue.get_resname() == "PTR":
                        residues_pdb.append("y")
                        residues_seq.append("Y")
                    # phosphorylated threonine
                    elif residue.get_resname() == "TPO":
                        residues_pdb.append("p")
                        residues_seq.append("T")
                    # phosphorylated serine
                    elif residue.get_resname() == "SEP":
                        residues_pdb.append("s")
                        residues_seq.append("S")
                    # seleno-methionine
                    elif residue.get_resname() == "MSE":
                        # no matter foldx, keep a methionine
                        residues_pdb.append("M")
                        residues_seq.append("M")
                    else:
                        residues_pdb.append(seq1(residue.get_resname()))
                        residues_seq.append(seq1(residue.get_resname()))
                    if len(residues_seq) == 1:
                        peptide_start = residue.id[1]
                        previous_aa = residue.id[1]
                        continue

                # info_gap is a dictionary with the gap coordinates
                if residue.id[1] != previous_aa + 1 and residue.id[0] in [" ", "H_PTR", "H_TPO", "H_SEP", "H_MSE"] and len(residue.get_resname().strip()) == 3:
                    info_gap = {}
                    dif = residue.id[1] - previous_aa
                    info_gap["Start"] = previous_aa
                    info_gap["End"] = residue.id[1]
                    info_gap["Difference"] = dif
                    info_gaps.append(info_gap)

                previous_aa = residue.id[1]

            # as motifs are generally present in small chains, we look for them
            if len(residues_seq) < 50:
                seq_chain = "".join(residues_seq)
                pdb_chain = "".join(residues_pdb)
                # look for matches in the sequence chain, where residues are as the standard 20,
                # and then we build a pdb_match object with it and keep it
                for hit in re.finditer(motif_re, str(seq_chain)):
                    pdb_match = pdb_match_object("alphafold", chain.id, peptide_start, hit, structure, info_gaps, pdb_chain)
                    return pdb_match, file[:-4]

def pdb_ddG_stability_matrix(motif_name, pdb_match, motif_label, tmp_path, repaired_pdbs_path, pdb_db_path, cwd_path):
    """
    This function takes as input a pdb_match, and will generate the stability ddG matrix
    In: motif name and re in point-based format, pdb_match, temporary files path,
    repaired pdbs path, pdb crystal structures path
    Out: stability ddG matrix
    """
    # check if we have already calculated the matrix for that pdb_match, if so, we load the file with the matrix, prevent losing time
    if os.path.exists(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix-{pdb_match.PDB}_{pdb_match.chain}.tsv"):
        stab_matrix = pd.read_csv(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix-{pdb_match.PDB}_{pdb_match.chain}.tsv", sep="\t", index_col="Motif")
        return stab_matrix

    else:
        # check if the crystal structure has already been repaired, prevent losing time
        if not os.path.exists(f"{repaired_pdbs_path}{pdb_match.PDB}_Repair.pdb"):
            # FoldX can only make calculations if the crystal structure file is in the working directory
            # If the file is not here, we look for it in the database and write a pdb file in the current working directory
            if not os.path.exists(f"{cwd_path}{pdb_match.PDB}.pdb"):
                archivo = gzip.open(f"{pdb_db_path}pdb{pdb_match.PDB.lower()}.ent.gz", "rb")
                pdb = archivo.read()
                with open(f"{cwd_path}{pdb_match.PDB}.pdb", "wb") as handle:
                    handle.write(pdb)
                    handle.close()
            print(f"Repairing {pdb_match.PDB} crystal")
            # repair the crystal structure, outputs all the files in tmp and then we move it to the working directory
            os.system(f"./foldx --command=RepairPDB --pdb={pdb_match.PDB}.pdb --output-dir={tmp_path} > {tmp_path}output_foldx.txt")
            os.system(f"mv {tmp_path}{pdb_match.PDB}_Repair.pdb {cwd_path}")
            os.system(f"rm {tmp_path}{pdb_match.PDB}_Repair.fxout")
            os.system(f"rm {pdb_match.PDB}.pdb")
        else:
            os.system(f"mv {repaired_pdbs_path}{pdb_match.PDB}_Repair.pdb {cwd_path}")

        # build positions in the crystal structure that we are going to mutate and analyze
        # we use motif_pdb sequence because is the residue code that foldx will always understand (think in phosphorylated residues)
        positions = []
        for i, aa in enumerate(pdb_match.motif_pdb):
            positions.append(f"{aa}{pdb_match.chain}{str(pdb_match.start_position+i)}a")
        positions = ",".join(positions)
        print(f"Running {pdb_match.PDB} PositionScan on chain {pdb_match.chain}")
        # calculate the stability ddG for each of the 19 alternative residues
        os.system(f"./foldx --command=PositionScan --pdb={pdb_match.PDB}_Repair.pdb --positions={positions} --output-dir={tmp_path} > {tmp_path}output_foldx.txt")
        # we dont need anymore the repaired pdb file, we move it with the others repaired pdbs
        if not os.path.exists(f"{repaired_pdbs_path}{pdb_match.PDB}_Repair.pdb"):
            os.system(f"mv {pdb_match.PDB}_Repair.pdb {repaired_pdbs_path}")
        # load the matrix as a pandas dataframe, the predictions for each residue are placed every 21 rows
        dfs = []
        n = 21
        # iterate through residues predictions
        for i in motif_label.split("."):
            df = pd.read_csv(f"{tmp_path}PS_{pdb_match.PDB}_Repair_scanning_output.txt", sep="\t", header=None)
            # we will take a multiple of 21 for each iteration corresponding to each residue
            df = df[n-21:n]
            df.reset_index(inplace=True, drop=True)
            # eliminate the first row (change by itself)
            df.drop(0, inplace=True)
            # transpose the matrix
            df = df.T
            # eliminate the first row with the changes, keeps the order every time (see foldx_aminoacids)
            df.drop(0, inplace=True)
            # keep the values as a single dataframe row
            dfs.append(df)
            # we sum 21 to capture the predictions of the next residue
            n += 21

        # concat all residues and set type of data
        stab_matrix = pd.concat(dfs)
        stab_matrix = stab_matrix[stab_matrix.columns].astype(float)
        # set the columns names as given by foldx
        foldx_aminoacids = ["G", "A", "L", "V", "I", "P", "R", "T", "S", "C", "M", "K", "E", "Q", "D", "N", "W", "Y", "F", "H"]
        stab_matrix.columns = foldx_aminoacids
        # change columns order in a more physiochemical way
        aminoacids = ["G", "A", "V", "L", "I", "M", "P", "W", "F", "Y", "S", "T", "C", "N", "Q", "D", "E", "K", "R", "H"]
        stab_matrix = stab_matrix[aminoacids]
        # round
        stab_matrix = round(stab_matrix, 1).replace(-0, 0)
        # index column with the motif
        stab_matrix["Motif"] = [f"{pos}_{n}" for n, pos in enumerate(motif_label.split("."))]
        stab_matrix.set_index("Motif", inplace=True)
        stab_matrix.to_csv(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix-{pdb_match.PDB}_{pdb_match.chain}.tsv", sep="\t")

        # remove some extra files
        aminoacids = ["GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PRO", "TRP", "PHE", "TYR", "SER", "THR", "CYS", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"]
        for i in range(len(motif_label.split("."))):
            for aa in aminoacids:
                os.system(f"rm {cwd_path}{aa}{pdb_match.start_position+i}_{pdb_match.PDB}_Repair.pdb")
            if os.path.exists(f"{cwd_path}H1S{pdb_match.start_position+i}_{pdb_match.PDB}_Repair.pdb"):
                os.system(f"rm {cwd_path}H1S{pdb_match.start_position+i}_{pdb_match.PDB}_Repair.pdb")
            if os.path.exists(f"{cwd_path}H2S{pdb_match.start_position+i}_{pdb_match.PDB}_Repair.pdb"):
                os.system(f"rm {cwd_path}H2S{pdb_match.start_position+i}_{pdb_match.PDB}_Repair.pdb")
            os.system(f"rm {cwd_path}energies_{pdb_match.start_position+i}_{pdb_match.PDB}_Repair.txt")
            if os.path.exists(f"{cwd_path}binding_energies_{pdb_match.start_position+i}_{pdb_match.PDB}_Repair.txt"):
                os.system(f"rm {cwd_path}binding_energies_{pdb_match.start_position+i}_{pdb_match.PDB}_Repair.txt")
        os.system(f"rm {tmp_path}PS_{pdb_match.PDB}_Repair_scanning_output.txt")
        os.system(f"rm {tmp_path}PS_{pdb_match.PDB}_Repair.fxout")

        return stab_matrix

def alphafold_ddG_stability_matrix(motif_name, pdb_match, row, file, motif_label, tmp_path, repaired_pdbs_path, cwd_path):
    """
    This function takes as input a pdb_match, and will generate the stability ddG matrix
    In: motif name and re in point-based format, pdb_match, temporary files path,
    repaired pdbs path, pdb crystal structures path
    Out: stability ddG matrix
    """
    # check if we have already calculated the matrix for that pdb_match, if so, we load the file with the matrix, prevent losing time
    ranking = row["rank"]
    if os.path.exists(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix_AF2_{motif_name}_rank_{ranking}_model_{row.model}.tsv"):
        stab_matrix = pd.read_csv(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix_AF2_{motif_name}_rank_{ranking}_model_{row.model}.tsv", sep="\t", index_col="Motif")
        return stab_matrix

    else:
        # check if the crystal structure has already been repaired, prevent losing time
        if not os.path.exists(f"{repaired_pdbs_path}{file}_Repair.pdb"):
            # FoldX can only make calculations if the crystal structure file is in the working directory
            # If the file is not here, we look for it in the database and write a pdb file in the current working directory
            if not os.path.exists(f"{cwd_path}{file}.pdb"):
                archivo = f"{cwd_path}Motifs/{motif_name}/Seed/{file}.pdb"
                os.system(f"cp {archivo} {cwd_path}")
            print(f"Repairing AF2 crystal")
            # repair the crystal structure, outputs all the files in tmp and then we move it to the working directory
            os.system(f"./foldx --command=RepairPDB --pdb={file}.pdb --output-dir={tmp_path} > {tmp_path}output_foldx.txt")
            os.system(f"mv {tmp_path}{file}_Repair.pdb {cwd_path}")
            os.system(f"rm {tmp_path}{file}_Repair.fxout")
            os.system(f"rm {cwd_path}{file}.pdb")
        else:
            os.system(f"mv {repaired_pdbs_path}{file}_Repair.pdb {cwd_path}")

        # build positions in the crystal structure that we are going to mutate and analyze
        # we use motif_pdb sequence because is the residue code that foldx will always understand (think in phosphorylated residues)
        positions = []
        for i, aa in enumerate(pdb_match.motif_pdb):
            positions.append(f"{aa}{pdb_match.chain}{str(pdb_match.start_position+i)}a")
        positions = ",".join(positions)
        print(file)
        print(positions)
        print(f"Running AF2 PositionScan on chain {pdb_match.chain}")
        # calculate the stability ddG for each of the 19 alternative residues
        os.system(f"./foldx --command=PositionScan --pdb={file}_Repair.pdb --positions={positions} --output-dir={tmp_path} > {tmp_path}output_foldx.txt")
        # we dont need anymore the repaired pdb file, we move it with the others repaired pdbs
        if not os.path.exists(f"{repaired_pdbs_path}{file}_Repair.pdb"):
            os.system(f"mv {file}_Repair.pdb {repaired_pdbs_path}")
        # load the matrix as a pandas dataframe, the predictions for each residue are placed every 21 rows
        dfs = []
        n = 21
        # iterate through residues predictions
        for i in motif_label.split("."):
            df = pd.read_csv(f"{tmp_path}PS_{file}_Repair_scanning_output.txt", sep="\t", header=None)
            # we will take a multiple of 21 for each iteration corresponding to each residue
            df = df[n-21:n]
            df.reset_index(inplace=True, drop=True)
            # eliminate the first row (change by itself)
            df.drop(0, inplace=True)
            # transpose the matrix
            df = df.T
            # eliminate the first row with the changes, keeps the order every time (see foldx_aminoacids)
            df.drop(0, inplace=True)
            # keep the values as a single dataframe row
            dfs.append(df)
            # we sum 21 to capture the predictions of the next residue
            n += 21

        # concat all residues and set type of data
        stab_matrix = pd.concat(dfs)
        stab_matrix = stab_matrix[stab_matrix.columns].astype(float)
        # set the columns names as given by foldx
        foldx_aminoacids = ["G", "A", "L", "V", "I", "P", "R", "T", "S", "C", "M", "K", "E", "Q", "D", "N", "W", "Y", "F", "H"]
        stab_matrix.columns = foldx_aminoacids
        # change columns order in a more physiochemical way
        aminoacids = ["G", "A", "V", "L", "I", "M", "P", "W", "F", "Y", "S", "T", "C", "N", "Q", "D", "E", "K", "R", "H"]
        stab_matrix = stab_matrix[aminoacids]
        # round
        stab_matrix = round(stab_matrix, 1).replace(-0, 0)
        # index column with the motif
        stab_matrix["Motif"] = [f"{pos}_{n}" for n, pos in enumerate(motif_label.split("."))]
        stab_matrix.set_index("Motif", inplace=True)
        stab_matrix.to_csv(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix_AF2_{motif_name}_rank_{ranking}_model_{row.model}.tsv", sep="\t")

        # remove some extra files
        aminoacids = ["GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PRO", "TRP", "PHE", "TYR", "SER", "THR", "CYS", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"]
        for i in range(len(motif_label.split("."))):
            for aa in aminoacids:
                os.system(f"rm {cwd_path}{aa}{pdb_match.start_position+i}_{file}_Repair.pdb")
            if os.path.exists(f"{cwd_path}H1S{pdb_match.start_position+i}_{file}_Repair.pdb"):
                os.system(f"rm {cwd_path}H1S{pdb_match.start_position+i}_{file}_Repair.pdb")
            if os.path.exists(f"{cwd_path}H2S{pdb_match.start_position+i}_{file}_Repair.pdb"):
                os.system(f"rm {cwd_path}H2S{pdb_match.start_position+i}_{file}_Repair.pdb")
            os.system(f"rm {cwd_path}energies_{pdb_match.start_position+i}_{file}_Repair.txt")
            if os.path.exists(f"{cwd_path}binding_energies_{pdb_match.start_position+i}_{file}_Repair.txt"):
                os.system(f"rm {cwd_path}binding_energies_{pdb_match.start_position+i}_{file}_Repair.txt")
        os.system(f"rm {tmp_path}PS_{file}_Repair_scanning_output.txt")
        os.system(f"rm {tmp_path}PS_{file}_Repair.fxout")

        return stab_matrix

def pdb_mean_foldx_matrix(total_pdb_matchs, motif_name, motif_label, cwd_path, plot=False):
    """
    Outputs the mean FoldX Matrix for a motif. It considers mainly the ddGbinding energies but if a significant
    destabilization effect is detected, it keeps the stability values for that residue
    In: list of pdb_matches where FoldX calculations were performed
    Out: mean ddG matrix and plot if wanted
    """
    # load all stability matrices and sum their values
    df = pd.read_csv(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix-{total_pdb_matchs[0].PDB}_{total_pdb_matchs[0].chain}.tsv", sep="\t", index_col=0)
    for pdb_match in total_pdb_matchs[1:]:
        df += pd.read_csv(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix-{pdb_match.PDB}_{pdb_match.chain}.tsv", sep="\t", index_col=0)
    # get the mean by dividing by the amount of pdb_matches
    stab_matrix = df/len(total_pdb_matchs)
    # round
    stab_matrix = round(stab_matrix, 1).replace(-0, 0)
    stab_matrix.reset_index(inplace=True)
    stab_matrix["Motif"] = motif_label.split(".")
    stab_matrix.set_index("Motif", inplace=True)
    stab_matrix.to_csv(f"{cwd_path}Motifs/{motif_name}/MeanPositionScanMatrix_pdb_{motif_name}.tsv", sep="\t")

    # plot the stability matrix
    sns.set(font_scale=1.2)
    plt.figure(figsize=(15, 5))
    cmap = sns.color_palette("PiYG_r", as_cmap=True)
    ax = sns.heatmap(stab_matrix, vmin=-5, vmax=5, cmap=cmap, annot=True, linecolor="white", linewidths=.5, cbar_kws={"label":"Kcal/mol"})
    plt.yticks(rotation=0)
    plt.ylabel("")
    plt.savefig(f"{cwd_path}Motifs/{motif_name}/MeanPositionScanMatrix_pdb_{motif_name}.png")
    if plot == True:
        plt.show()
    else:
        plt.close()

    return stab_matrix

def alphafold_mean_foldx_matrix(total_pdb_matchs, filtered_predictions_df, motif_name, motif_label, cwd_path, plot=False):
    """
    Outputs the mean FoldX Matrix for a motif. It considers mainly the ddGbinding energies but if a significant
    destabilization effect is detected, it keeps the stability values for that residue
    In: list of pdb_matches where FoldX calculations were performed
    Out: mean ddG matrix and plot if wanted
    """

    # load all stability matrices and sum their values
    switch = False
    for n, row in filtered_predictions_df.iterrows():
        ranking = row["rank"]
        if switch == False:
            df = pd.read_csv(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix_AF2_{motif_name}_rank_{ranking}_model_{row.model}.tsv", sep="\t", index_col=0)
            switch = True
        elif switch == True:
            df += pd.read_csv(f"{cwd_path}Motifs/{motif_name}/PositionScanMatrix_AF2_{motif_name}_rank_{ranking}_model_{row.model}.tsv", sep="\t", index_col=0)

    # get the mean by dividing by the amount of pdb_matches
    stab_matrix = df/len(total_pdb_matchs)
    # round
    stab_matrix = round(stab_matrix, 1).replace(-0, 0)
    stab_matrix.reset_index(inplace=True)
    stab_matrix["Motif"] = motif_label.split(".")
    stab_matrix.set_index("Motif", inplace=True)
    stab_matrix.to_csv(f"{cwd_path}Motifs/{motif_name}/MeanPositionScanMatrix_AF2_{motif_name}.tsv", sep="\t")

    # plot the stability matrix
    sns.set(font_scale=1.2)
    plt.figure(figsize=(15, 5))
    cmap = sns.color_palette("PiYG_r", as_cmap=True)
    ax = sns.heatmap(stab_matrix, vmin=-5, vmax=5, cmap=cmap, annot=True, linecolor="white", linewidths=.5, cbar_kws={"label":"Kcal/mol"})
    plt.yticks(rotation=0)
    plt.ylabel("")
    plt.savefig(f"{cwd_path}Motifs/{motif_name}/MeanPositionScanMatrix_AF2_{motif_name}.png")
    if plot == True:
        plt.show()
    else:
        plt.close()

    return stab_matrix

def foldx_process(true_positives, motif_re, motif_name, motif_label, tmp_path, repaired_pdbs_path, pdb_db_path, cwd_path, plot = False):
    """
    This function drives the general process of building the PositionScan substitution matrices
    """

    print("Building Motif's PositionScan Substitution Matrices...")

    # list where storing pdb_match objects
    total_pdb_matchs = []

    if len([match for match in true_positives if len(match.PDBs) > 0]) > 0:
        source = "pdb"
        # iterate through true positives motif objects that contain crystal structures on them and get the
        # pdb_match objects that define interactions
        for match in tqdm([match for match in true_positives if len(match.PDBs) > 0]):
            pdb_matchs = pdb_foldx_parameters(match, motif_re, motif_name, pdb_db_path)
            # iterate through pdb_matches and calculate ddGb with FoldX while appending them to the final list
            for pdb_match in pdb_matchs:
                #ddGb_matrix = pdb_ddG_binding_matrix(motif_name, pdb_match, motif_label, tmp_path, repaired_pdbs_path, pdb_db_path, cwd_path)
                stab_matrix = pdb_ddG_stability_matrix(motif_name, pdb_match, motif_label, tmp_path, repaired_pdbs_path, pdb_db_path, cwd_path)
                total_pdb_matchs.append(pdb_match)

        # calculate the final mean matrix
        MeanFoldXMatrix = pdb_mean_foldx_matrix(total_pdb_matchs, motif_name, motif_label, cwd_path, plot=plot)

    else:
        source = "af2"
        predictions_df = pd.read_csv(f"{cwd_path}Motifs/{motif_name}/interface_res_stat_total.csv")
        filtered_predictions_df = predictions_df[(predictions_df.Interacion_energy <= -5) & (predictions_df.DDT_mean_motif_IR >= 65)]
        if filtered_predictions_df.empty:
            print("Reliable stability matrices could not be generated")
            return None, source
        for n, row in filtered_predictions_df.iterrows():
            pdb_match, file = alphafold_foldx_parameters(motif_re, motif_name, row, cwd_path)
            #ddGb_matrix = alphafold_ddG_binding_matrix(motif_name, pdb_match, motif_label, tmp_path, repaired_pdbs_path, pdb_db_path, cwd_path)
            stab_matrix = alphafold_ddG_stability_matrix(motif_name, pdb_match, row, file, motif_label, tmp_path, repaired_pdbs_path, cwd_path)
            total_pdb_matchs.append(pdb_match)

        MeanFoldXMatrix = alphafold_mean_foldx_matrix(total_pdb_matchs, filtered_predictions_df, motif_name, motif_label, cwd_path, plot=plot)

    return MeanFoldXMatrix, source

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Returns a list of match objects for the motif on a sequences database (default SwissProt)")
    parser.add_argument("Motif_name", action="store", help = "Name of the motif, e.g. DOC_MAPK_JIP1_4")
    parser.add_argument("Motif_re", action="store", help = "Motif regular expression, e.g. [RK]P[^P][^P]L.[LIVMF]")
    parser.add_argument("Motif_label", action="store", help = "Motif regular expression for plotting, e.g. RK.P.^P.^P.L.x.LIVMF")
    parser.add_argument("-fasta", "--sp_h_fasta", action="store", help = "file with fasta entries for SwissProt", default="./UniProt/uniprot_sprot_h.fasta")
    parser.add_argument("-xml", "--sp_h_xml", action="store", help = "file with xml entries for SwissProt", default="./UniProt/uniprot_sprot_h.xml")
    parser.add_argument("-list", "--sp_h_list", action="store", help = "file with all UniProtIDs in SwissProt", default="./UniProt/uniprot_sprot_h.list")
    parser.add_argument("-elm_inst", "--elm_instances", action="store", help = "ELM Instances database where look for motifs", default="./ELM/elm_instances.tsv")
    parser.add_argument("-pdb_str", "--pdb_structures", action="store", help = "PDB database with all crystals", default="./PDB/zipped/")
    parser.add_argument("-rep_pdbs", "--repaired_pdbs", action="store", help = "folder with repaired PDBSs", default="./repaired_pdbs/")
    parser.add_argument("-tmp", "--tmp_path", action="store", help = "temporary folder", default="./tmp/")
    parser.add_argument("-cwd", "--cwd", action="store", help = "current working directory path", default="./")

    args = parser.parse_args()
    
    if not os.path.exists(f"{args.cwd}Motifs/{args.Motif_name}"):
        os.system(f"mkdir {args.cwd}Motifs/{args.Motif_name}")

    if not os.path.exists(f"{args.tmp_path}"):
        os.system(f"mkdir {args.tmp_path}")

    from MotifSearcher import HumanProteomeSearch
    screening_proteoma = HumanProteomeSearch(args.Motif_name, args.Motif_re, args.sp_h_fasta, args.sp_h_xml, args.sp_h_list)

    from positive_ctrl import positive_ctrl
    positive_ctrl(screening_proteoma, args.Motif_name, args.pdb_structures, args.elm_instances)

    true_positives = [match for match in screening_proteoma if match.ELM_Positive == True]
    FoldXMatrix = foldx_process(true_positives, args.Motif_re, args.Motif_name, args.Motif_label, args.tmp_path, args.repaired_pdbs, args.pdb_structures, args.cwd)
