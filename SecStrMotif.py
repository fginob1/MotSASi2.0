#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
from Bio import SeqIO
import os
import pandas as pd
import time
import ssbio.protein.sequence.properties.scratch
import numpy as np
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import argparse
from Bio.PDB import PDBIO, PDBParser
import gzip
from Bio.PDB.PDBIO import Select
import biotite.structure.io as strucio
import biotite.structure as struc

def prot_fasta(UniProtID, seq, sec_str_path):
    """
    Esta funcion toma el UniProtID de la proteina en estudio y su correspondiente objeto secuencia
    y arma un archivo fasta que se utilizara posteriormente para realizar blast contra pdb
    """
    with open(f"{sec_str_path}tmp_file_{UniProtID}.fasta", "w") as out_handle:
        out_handle.write(f">{UniProtID}\n{str(seq)}")
        out_handle.close()
    return

class Chain_Select(Select):
    """
    Auxiliar function to isolate a specific chain in a new pdb file
    """
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, cadena):
        return self.chain == cadena.id

def pred_est_sec_x_residuo(ctrl_positivo, fasta_tmp_file, seq, UniProtID, sec_str_path, alphafold_human_path, scratch_path, tmp_path, pred_number):
    """
    This function takes as input a gene and its protein product (as given by UniProt) and outputs its
    secondary structure prediction
    """
    
    # we first try to take secondary structure from alphafold
    p = PDBParser()
    structure_id = UniProtID
    structure = p.get_structure(structure_id, gzip.open(f"{alphafold_human_path}AF-{ctrl_positivo.UniProtID}-F1-model_v4.pdb.gz", "rt"))
    chain = "A"
    io = PDBIO()
    io.set_structure(structure)
    io.save(f"{tmp_path}chain_only_aux_{ctrl_positivo.UniProtID}.pdb", select = Chain_Select(chain = chain))
    p = PDBParser()
    filename = f"{tmp_path}chain_only_aux_{ctrl_positivo.UniProtID}.pdb"
    structure = p.get_structure(structure_id, filename)
    array = strucio.load_structure(filename)
    sse = struc.annotate_sse(array)

    # lists with dssp numerations, residues and secondary structure motif secstr_motifs, respectively
    residue_number = [n+1 for n, aa in enumerate(seq)]
    residues = [aa for aa in seq]
    verdicts = list(sse)
    verdicts = list(map(lambda x: x.replace("a", "H"), verdicts))
    verdicts = list(map(lambda x: x.replace("b", "E"), verdicts))
    verdicts = list(map(lambda x: x.replace("c", "C"), verdicts))
    confidences = [atom.bfactor for model in structure for chain in model for residue in chain for atom in residue if atom.name == "CA"]

    check = False
    # # check if the length of the sequences of UniProt and AlphaFold are the same
    if len(residues) == len(verdicts):
        # dictionary with lists and pandas dataframe with all the info
        secstr_raw = {"Residue_number": residue_number, "Residue": residues, "SecStr": verdicts, "Confidence": confidences}
        secstr_df = pd.DataFrame(secstr_raw, columns = ["Residue_number", "Residue", "SecStr", "Confidence"])
        secstr_df.set_index("Residue_number", inplace=True)
        # check if the length of the sequences are the same and if the alphafold confidence is good
        #if len(seq) == len(residues) and secstr_df[ctrl_positivo.start:ctrl_positivo.start+len(ctrl_positivo.seq)]["Confidence"].min() > 80:
        if secstr_df[ctrl_positivo.start:ctrl_positivo.start+len(ctrl_positivo.seq)]["Confidence"].min() > 80:
            check = True

    if check == True:
        return secstr_df
    else:
        # SCRATCH class
        resultado = ssbio.protein.sequence.properties.scratch.SCRATCH(UniProtID, seq_file = fasta_tmp_file)
        # method to generate the prediction, if already done parse it automatically
        trabajo = resultado.run_scratch(f"{scratch_path}run_SCRATCH-1D_predictors.sh", num_cores = 4, outname=UniProtID, outdir=f"{sec_str_path}", force_rerun=False)

        # method to parse the result, it is a dictionary with the UniProtID as key and the prediction as value
        if pred_number == 3:
            prediccion = resultado.sspro_results()[UniProtID]
        elif pred_number == 8:
            prediccion = resultado.sspro8_results()[UniProtID]

        # list with prediction by residue
        pred_est_sec = [i for i in prediccion]

        # lists with residues and their numeration (1-based)
        residues = []
        residue_number = []
        for c, i in enumerate(seq):
            residues.append(i)
            residue_number.append(c+1)

        # dictionary with lists as values
        Residuos_est_sec= {"Residue_number": residue_number, "Residue": residues, "SecStr": pred_est_sec}

        # pandas dataframe from the dictionary and set residue number as index
        df_est_sec = pd.DataFrame(Residuos_est_sec, columns = ["Residue_number", "Residue", "SecStr"])
        df_est_sec.set_index("Residue_number", inplace=True)

        return df_est_sec

def scratch(matchs_list, motif_name, motif_label, sp_h_xml_path, sec_str_path, alphafold_human_path, scratch_path, cwd_path, tmp_path, plot = True, screen = False, save_tsv = True, pred_number = 3):
    """
    This function takes as input the positive control group (or a specific protein) and outputs a dataframe
    describing the percentage of secondary structure presentations for each residue in the motif
    """

    # parse UniProt, SwissProt database in xml format
    xmls = SeqIO.index_db(f"{sp_h_xml_path}.idx", sp_h_xml_path, "uniprot-xml")

    # length of the motif
    M_length = len(motif_label.split("."))

    # dictionary to keep secondary structure info, keys are UniProtIDs followed by start position of the motif,
    # values are a list of length equal to the motif length with residues secondary structure predictions
    diccionario_est_sec = {}

    if len(matchs_list) == 1:
        # make the same iteration but without tqdm, same code
        for ctrl_positivo in matchs_list:
            # parse UniProt, SwissProt entry
            seq_record = xmls[ctrl_positivo.UniProtID]
            # UniProt protein sequence
            seq = seq_record.seq
            # generate the fasta file with the protein sequence
            prot_fasta(ctrl_positivo.UniProtID, seq, sec_str_path)
            # make the secondary structure prediction and generate a pandas dataframe with protein residues
            # and their corresponding secondary structure
            prediccion = pred_est_sec_x_residuo(ctrl_positivo, f"{sec_str_path}tmp_file_{ctrl_positivo.UniProtID}.fasta", seq, ctrl_positivo.UniProtID, sec_str_path, alphafold_human_path, scratch_path, tmp_path, pred_number)
            # generate a list with np.nan values that we are going to replace by predictions white iterating
            lista_aux = [np.nan]*M_length
            for n in range(len(ctrl_positivo.seq)):
                lista_aux[n] = prediccion.loc[ctrl_positivo.start + n].SecStr
            # keep the list in the dictionary with UniProtID plus start position as key
            diccionario_est_sec[ctrl_positivo.UniProtID + "_" + str(ctrl_positivo.start)] = lista_aux
            os.remove(f"{sec_str_path}tmp_file_{ctrl_positivo.UniProtID}.fasta")
    else:
        # iterate through motifs in the positive control group
        for ctrl_positivo in tqdm(matchs_list):
            # parse UniProt, SwissProt entry
            seq_record = xmls[ctrl_positivo.UniProtID]
            # UniProt protein sequence
            seq = seq_record.seq
            # generate the fasta file with the protein sequence
            prot_fasta(ctrl_positivo.UniProtID, seq, sec_str_path)
            # make the secondary structure prediction and generate a pandas dataframe with protein residues
            # and their corresponding secondary structure
            prediccion = pred_est_sec_x_residuo(ctrl_positivo, f"{sec_str_path}tmp_file_{ctrl_positivo.UniProtID}.fasta", seq, ctrl_positivo.UniProtID, sec_str_path, alphafold_human_path, scratch_path, tmp_path, pred_number)
            # generate a list with np.nan values that we are going to replace by predictions white iterating
            lista_aux = [np.nan]*M_length
            for n in range(len(ctrl_positivo.seq)):
                lista_aux[n] = prediccion.loc[ctrl_positivo.start + n].SecStr
            # keep the list in the dictionary with UniProtID plus start position as key
            diccionario_est_sec[ctrl_positivo.UniProtID + "_" + str(ctrl_positivo.start)] = lista_aux
            os.remove(f"{sec_str_path}tmp_file_{ctrl_positivo.UniProtID}.fasta")

    # generate a pandas dataframe from the dictionary, rows are motif residues, columns are proteins with predictions
    df = pd.DataFrame(diccionario_est_sec)
    # transpose the dataframe
    df = df.T
    # count the ocurrences of secondary structure predictions by column, giving as ouput a dataframe with
    # motif residues in columns and secondary structure motifs in rows
    df = df.apply(pd.Series.value_counts)
    # transpose, we obtain motif residues in rows and secondary structure motifs in columns
    df = df.T
    # replace empty cells by zeros
    df = df.fillna(0)
    # dictionary with the letter-secondary structure motif code for scratch (our predictor)
    sec_est = {"H": "alpha-helix", "G": "310-helix", "I": "pi-helix (extremely rare)", 
           "E": "extended strand", "B": "beta-bridge", "T": "turn", "S": "bend", "C": "the rest"}
    # replace explicit motifs by letters
    df = df.rename(columns=lambda x: sec_est[x])
    # sum total events by row
    df_total = df.sum(axis=1, numeric_only=True)
    # convert to percentage of events
    df_rel = df[df.columns[:]].div(df_total, 0)*100
    # rows labels, create a column and set them as index
    motif_label = motif_label.split(".")
    df_rel["Motif"] = motif_label
    df_rel.set_index("Motif", inplace=True, drop=True)

    if save_tsv == True:
        df_rel.to_csv(f"{cwd_path}Motifs/{motif_name}/SecEst_{motif_name}.tsv", sep="\t")

    if plot == True:
        # invert the dataframe to get the residues in order in the plot
        df_rel = df_rel.iloc[::-1]
        df_rel.reset_index(inplace=True)
        # plotting
        df_rel.plot(x = "Motif", kind = "barh", stacked = True, figsize=(7,9), fontsize=14)
        plt.ylabel("Motif", fontsize=18)
        plt.title("Secondary Structure prediction", fontsize=18)
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Secondary Structure elements")
        # percentage labels in the plot
        # iterate through secondary structure motifs in columns (omit first column, motif residues)
        for n in df_rel.columns[1:]:
            # iterate through rows from which a zip was made between sec str motif in a residue and cumulated
            # sum of the percentages in that row
            for i, (pc, cs) in enumerate(zip(df_rel[n], df_rel.iloc[:, 1:].cumsum(1)[n])):
                # only incorporate label if ocurrence percentage is more than 5%
                if pc >= 5:
                    # x coordinate is equal to cumulated percentage minus half of the percentage of the sec str
                    # motif being evaluated, y coordinate is given by the residue evaluated
                    plt.text(cs-pc/2, i, f"{str(np.round(pc, 1))}%", va = "center", ha = "center", fontsize=14)

        plt.savefig(f"{cwd_path}Motifs/{motif_name}/SecEst_{motif_name}.png")
        if screen == True:
            plt.show()
        else:
            plt.close()

        # regenerate the index
        df_rel = df_rel.iloc[::-1]
        df_rel["Motif"] = motif_label
        df_rel.set_index("Motif", inplace=True, drop=True)

    return df_rel

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
    parser.add_argument("-sec_str", "--secondary_structures", action="store", help = "folder with secondary structure predictions", default="./secondary_structures/")
    parser.add_argument("-af_human", "--alphafold_human", action="store", help = "alphafold for the human proteome", default="./AlphaFold2_human_proteome/UP000005640_9606_HUMAN_v4/")
    parser.add_argument("-scratch", "--scratch_path", action="store", help = "scratch secondary structure predictor path", default="./scratch/SCRATCH-1D_1.3/bin/")
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
    print("Scanning Motif Secondary Structure in positive controls")
    SecStr_df = scratch(true_positives, args.Motif_name, args.Motif_label, args.sp_h_xml, args.secondary_structures, args.alphafold_human, args.scratch_path, args.cwd, args.tmp_path)
