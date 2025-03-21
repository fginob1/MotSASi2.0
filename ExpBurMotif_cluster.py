#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
from Bio import ExPASy
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
import freesasa
from Bio.PDB.Polypeptide import three_to_one
from Bio import Align
from Bio.Align import substitution_matrices

class Chain_Select(Select):
    """
    Auxiliar function to isolate a specific chain in a new pdb file
    """
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, cadena):
        return self.chain == cadena.id

def get_fragment(query_pos, seq_length, fragments_length=1400, overlap=200):
    fragments = []
    start = 1
    end = 0
    
    # Crear los fragments hasta cubrir todos los residuos de la proteína
    while end < seq_length:
        end = start + fragments_length - 1
        if end > seq_length:
            end = seq_length
        fragments.append((start, end))
        start += overlap

    # Encontrar el fragmento que contiene la posición de interés
    chosen_fragment = None
    minimal_central_distance = float('inf')

    for i, (start, end) in enumerate(fragments):
        # Si la posición de interés está en este fragmento
        if start <= query_pos <= end:
            fragment_center = (start + end) // 2
            central_distance = abs(fragment_center - query_pos)

            # Buscar el fragmento donde la posición esté más cerca del centro
            if central_distance < minimal_central_distance:
                minimal_central_distance = central_distance
                chosen_fragment = f"F{i+1}"
                fragment = (start, end)
    
    return chosen_fragment, fragment

def expent_x_residuo(UniProtID, start, sp_h_xml_path, sasa_thr_path, alphafold_human_path, tmp_path):
    """
    Este script toma toda la info generada previamente y va a generar un dataframe con los residuos de la
    proteina en las filas y para cada uno de ellos que caiga en un dominio pfam analiza dicho residuo en una
    cantidad de cristales que presentan homologia y emite un veredicto sobre su exposicion
    """

    # levanto la base de datos xml de UniProt, SwissProt
    uniprot = SeqIO.index_db(f"{sp_h_xml_path}.idx", sp_h_xml_path, "uniprot-xml")

    seq_record = uniprot[UniProtID]
    seq = seq_record.seq
    seq_length = len(seq)

    # parse structure from alphafold, considerations should be made based on the length of the protein
    p = PDBParser()
    structure_id = UniProtID
    if len(seq) <= 2700:
        structure = p.get_structure(structure_id, gzip.open(f"{alphafold_human_path}AF-{UniProtID}-F1-model_v4.pdb.gz", "rt"))
    else:
        fragment_name, fragment = get_fragment(start, seq_length)
        structure = p.get_structure(structure_id, gzip.open(f"{alphafold_human_path}AF-{UniProtID}-{fragment_name}-model_v4.pdb.gz", "rt"))
    chain_ = "A"
    io = PDBIO()
    io.set_structure(structure)
    io.save(f"{tmp_path}chain_only_aux_{UniProtID}_{start}.pdb", select = Chain_Select(chain = chain_))

    residues = []
    for model in structure:
        for chain in model:
            if chain.id == chain_:
                for residue in chain:
                    if residue.get_id()[0] == " ":
                        residues.append(three_to_one(residue.resname))

    # voy a meter todos los polipeptidos de esa cadena del cristal en una unica secuencia que despues
    # el alineador me los "mapea" contra el fasta
    raw_seq_crystal = "".join(residues)

    #alineamos
    seq1 = seq
    seq2 = raw_seq_crystal
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.query_left_open_gap_score = 0
    aligner.query_left_extend_gap_score = 0
    aligner.query_right_open_gap_score = 0
    aligner.query_right_extend_gap_score = 0
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(seq1, seq2)
    alignment = alignments[0]
    seq_UniProt_aln = alignment[0]
    seq_alphafold_aln = alignment[1]

    # load csv file with exposed/buried thresholds
    cut_offs = pd.read_csv(sasa_thr_path)

    structure = freesasa.Structure(f"{tmp_path}chain_only_aux_{UniProtID}_{start}.pdb")
    result = freesasa.calc(structure)
    residueAreas = result.residueAreas()

    verdicts = []

    alphafold_residue_number = 0
    UniProt_residue_number = 0
    # starts iterating through alphafold alignment sequence (both residues and gaps)
    for i, UniProt_residue in enumerate(seq_UniProt_aln):
        if UniProt_residue == "-":
            alphafold_residue_number += 1
        else:
            UniProt_residue_number += 1
            if seq_alphafold_aln[i] == "-":
                verdicts.append("Unknown")
            else:
                alphafold_residue_number += 1
                res_type = residueAreas[chain_][str(alphafold_residue_number)].residueType
                SASA = residueAreas[chain_][str(alphafold_residue_number)].relativeSideChain
                if SASA >= cut_offs[cut_offs["Residuo"] == res_type]["Corte"].tolist()[0]:
                    verdicts.append("Exposed")
                elif SASA < cut_offs[cut_offs["Residuo"] == res_type]["Corte"].tolist()[0]:
                    verdicts.append("Buried")
                else:
                    verdicts.append("Unknown")

    # lists with residues and their numeration (1-based)
    residues = []
    residue_number = []
    for c, i in enumerate(seq):
        residues.append(i)
        residue_number.append(c+1)
    
    residues_expent = {"Residue_number": residue_number, "Residue": residues, "ExpEnt": verdicts}

    exp_ent_df = pd.DataFrame(residues_expent, columns = ["Residue_number", "Residue", "ExpEnt"])
    exp_ent_df.set_index("Residue_number", inplace=True)

    return exp_ent_df


def sasa(matchs_list, motif_name, motif_label, sp_h_xml_path, alphafold_human_path, sasa_thr_path, cwd_path, tmp_path, plot = True, screen = False, save_tsv = True, pred_number = 3):
    """
    This function takes as input the positive control group (or a specific protein) and outputs a dataframe
    describing the percentage of secondary structure presentations for each residue in the motif
    """

    # parse UniProt, SwissProt database in xml format
    xmls = SeqIO.index_db(f"{sp_h_xml_path}.idx", sp_h_xml_path, "uniprot-xml")

    # length of the motif
    M_length = len(motif_label.split("."))

    # dictionary to keep exposure info, keys are UniProtIDs followed by start position of the motif,
    # values are a list of length equal to the motif length with residues exposure predictions
    exp_bur_dictionary = {}

    if len(matchs_list) == 1:
        # make the same iteration but without tqdm, same code
        for ctrl_positivo in matchs_list:
            # parse UniProt, SwissProt entry
            seq_record = xmls[ctrl_positivo.UniProtID]
            # UniProt protein sequence
            seq = seq_record.seq
            # make the secondary structure prediction and generate a pandas dataframe with protein residues
            # and their corresponding secondary structure
            prediccion = expent_x_residuo(ctrl_positivo.UniProtID, ctrl_positivo.start, sp_h_xml_path, sasa_thr_path, alphafold_human_path, tmp_path)
            # generate a list with np.nan values that we are going to replace by predictions white iterating
            lista_aux = [np.nan]*M_length
            for n in range(len(ctrl_positivo.seq)):
                lista_aux[n] = prediccion.loc[ctrl_positivo.start + n].ExpEnt
            # keep the list in the dictionary with UniProtID plus start position as key
            exp_bur_dictionary[ctrl_positivo.UniProtID + "_" + str(ctrl_positivo.start)] = lista_aux
    else:
        # iterate through motifs in the positive control group
        for ctrl_positivo in tqdm(matchs_list):
            # parse UniProt, SwissProt entry
            seq_record = xmls[ctrl_positivo.UniProtID]
            # UniProt protein sequence
            seq = seq_record.seq
            # make the secondary structure prediction and generate a pandas dataframe with protein residues
            # and their corresponding secondary structure
            prediccion = expent_x_residuo(ctrl_positivo.UniProtID, ctrl_positivo.start, sp_h_xml_path, sasa_thr_path, alphafold_human_path, tmp_path)
            # generate a list with np.nan values that we are going to replace by predictions white iterating
            lista_aux = [np.nan]*M_length
            for n in range(len(ctrl_positivo.seq)):
                lista_aux[n] = prediccion.loc[ctrl_positivo.start + n].ExpEnt
            # keep the list in the dictionary with UniProtID plus start position as key
            exp_bur_dictionary[ctrl_positivo.UniProtID + "_" + str(ctrl_positivo.start)] = lista_aux

    # generate a pandas dataframe from the dictionary, rows are motif residues, columns are proteins with predictions
    df = pd.DataFrame(exp_bur_dictionary)
    # transpose the dataframe
    df = df.T
    # count the ocurrences of exposure predictions by column, giving as ouput a dataframe with
    # motif residues in columns and secondary structure motifs in rows
    df = df.apply(pd.Series.value_counts)
    # transpose, we obtain motif residues in rows and exposure in columns
    df = df.T
    # replace empty cells by zeros
    df = df.fillna(0)
    # sum total events by row
    df_total = df.sum(axis=1, numeric_only=True)
    # convert to percentage of events
    df_rel = df[df.columns[:]].div(df_total, 0)*100
    # rows labels, create a column and set them as index
    motif_label = motif_label.split(".")
    df_rel["Motif"] = motif_label
    df_rel.set_index("Motif", inplace=True, drop=True)

    if save_tsv == True:
        df_rel.to_csv(f"{cwd_path}Motifs/{motif_name}/ExpBur_{motif_name}.tsv", sep="\t")

    if plot == True:
        # invert the dataframe to get the residues in order in the plot
        df_rel = df_rel.iloc[::-1]
        df_rel.reset_index(inplace=True)
        # plotting
        df_rel.plot(x = "Motif", kind = "barh", stacked = True, figsize=(7,9), fontsize=14)
        plt.ylabel("Motif", fontsize=18)
        plt.title("Exposure prediction", fontsize=18)
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Exposure code")
        # percentage labels in the plot
        # iterate through secondary structure motifs in columns (omit first column, motif residues)
        for n in df_rel.columns[1:]:
            # iterate through rows from which a zip was made between sec str motif in a residue and cumulated
            # sum of the percentages in that row
            for i, (pc, cs) in enumerate(zip(df_rel[n], df_rel.iloc[:, 1:].cumsum(1)[n])):
                # only incorporate label if ocurrence percentage is more than 5%
                if pc >= 5:
                    # x coordinate is equal to cumulated percentage minus half of the percentage of the exposure
                    # being evaluated, y coordinate is given by the residue evaluated
                    plt.text(cs-pc/2, i, f"{str(np.round(pc, 1))}%", va = "center", ha = "center", fontsize=14)

        plt.savefig(f"{cwd_path}Motifs/{motif_name}/ExpBur_{motif_name}.png")
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
    parser.add_argument("-fasta", "--sp_h_fasta", action="store", help = "file with fasta entries for SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.fasta")
    parser.add_argument("-xml", "--sp_h_xml", action="store", help = "file with xml entries for SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.xml")
    parser.add_argument("-list", "--sp_h_list", action="store", help = "file with all UniProtIDs in SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.list")
    parser.add_argument("-elm_inst", "--elm_instances", action="store", help = "ELM Instances database where look for motifs", default="/grupos/Marce/estructural/databases/ELM/elm_instances.tsv")
    parser.add_argument("-pdb_str", "--pdb_structures", action="store", help = "PDB database with all crystals", default="/grupos/Marce/estructural/databases/PDB/zipped/")
    parser.add_argument("-exp_thr", "--exposure_thresholds", action="store", help = "exposure cut-offs", default="/grupos/Marce/estructural/tools/SASA/Cortes_SASA.csv")
    parser.add_argument("-af_human", "--alphafold_human", action="store", help = "alphafold for the human proteome", default="/grupos/Marce/estructural/databases/AlphaFold2/UP000005640_9606_HUMAN_v4/")
    parser.add_argument("-tmp", "--tmp_path", action="store", help = "temporary folder", default="/grupos/Marce/estructural/motsasi/tmp/")
    parser.add_argument("-cwd", "--cwd", action="store", help = "current working directory path", default="/grupos/Marce/estructural/motsasi/")

    args = parser.parse_args()
    
    if not os.path.exists(f"{args.cwd}Motifs/{args.Motif_name}"):
        os.system(f"mkdir {args.cwd}Motifs/{args.Motif_name}")

    if not os.path.exists(f"{args.tmp_path}"):
        os.system(f"mkdir {args.tmp_path}")

    from MotifSearcher_cluster import HumanProteomeSearch
    screening_proteoma = HumanProteomeSearch(args.Motif_name, args.Motif_re, args.sp_h_fasta, args.sp_h_xml, args.sp_h_list)

    from positive_ctrl_cluster import positive_ctrl
    positive_ctrl(screening_proteoma, args.Motif_name, args.pdb_structures, args.elm_instances)

    true_positives = [match for match in screening_proteoma if match.ELM_Positive == True]
    print("Scanning Motif Exposure in positive controls")
    sasa(true_positives, args.Motif_name, args.Motif_label, args.sp_h_xml, args.alphafold_human, args.exposure_thresholds, args.cwd, args.tmp_path)
