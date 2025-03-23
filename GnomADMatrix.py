#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
import math
from tqdm import tqdm
import warnings
import os
import argparse

warnings.filterwarnings("ignore")

def MatSubGnomAD(matchs_list, M_length, aas):
    """
    Compiles high allele frequency changes of each amino acid in the motif
    """

    # levanto el archivo con las variantes missense en GnomAD y su frecuencia
    #var_GnomAD = pd.read_csv(gnomad_missense_path)
    
    # Dictionary of dictionaries
    # We initialize a dictionary where keys are the 20 standard residues in the order given
    # Every residue of these ones has as value a dictionary where keys are indices of the positions
    # in the regular expression that also have an order, which is important in order to store the information
    # in an ordered way in lists
    diccionario_final = {}
    for aa in aas:
        diccionario_final[aa] = {}
        for i in range(M_length):
            diccionario_final[aa][i] = []

    # iterate through positive control matches
    for ctrl_positivo in tqdm(matchs_list):
        # iterate through previously generated variants
        for variante_GnomAD in [variant for variant in ctrl_positivo.Variants if variant.source_db == "GnomAD" and variant.verdict == "Benign"]:
            #  we enter into the dictionary the AF following the regular expression index and the residue change
            diccionario_final[variante_GnomAD.alt][variante_GnomAD.index].append(variante_GnomAD.confidence)

    return diccionario_final


def get_GnomADMatrix(matchs_list, motif_label, motif_name, cwd_path, iteration=False, number=False, plot = False):
    """
    This function takes as  input the motifs included in the positive control group and generates a matrix
    showing the reported aminoacidic changes with high allele frequency (AF) in GnomAD
    We can run this script both with pipeline_1 as also while running the iteration process
    """

    print("Building Motif's AF Substitution Matrix...")

    # list with aminoacids ordered in a physiochemical way
    aas = ["G", "A", "V", "L", "I", "M", "P", "W", "F", "Y", "S", "T", "C", "N", "Q", "D", "E", "K", "R", "H"]
    
    # length of the motif
    M_length = len(motif_label.split("."))
    
    # call the function that generates the dictionary that collects all the GnomAD data
    diccionario_final = MatSubGnomAD(matchs_list, M_length, aas)

    # generate a pandas dataframe with the dictionary, residues of the motif in rows,
    # 20 standard residues in columns
    GnomAD_matrix = pd.DataFrame(diccionario_final, columns=aas)
    
    # rows labels, create a column and set them as index
    motif_label = motif_label.split(".")
    GnomAD_matrix["Residuos"] = motif_label
    GnomAD_matrix.set_index("Residuos", inplace=True)
    
    # name the columns and save (we have lists in cells)
    GnomAD_matrix.columns = aas
    GnomAD_matrix.to_csv(f"{cwd_path}Motifs/{motif_name}/GnomADMatrix_{motif_name}_positives.tsv", sep="\t")

    # calculate mean values by list
    for aa in diccionario_final:
        for n in range(M_length):
            if len(diccionario_final[aa][n]) > 0:
                diccionario_final[aa][n] = max(diccionario_final[aa][n], key=abs)
            else:
                diccionario_final[aa][n] = np.nan

    # generate a pandas dataframe
    mean_GnomAD_matrix = pd.DataFrame(diccionario_final, columns=aas)

    # rows labels, create a column and set them as index
    mean_GnomAD_matrix["Residuos"] = motif_label
    mean_GnomAD_matrix.set_index("Residuos", inplace=True)

    # name columns
    mean_GnomAD_matrix.columns = aas

    # save the dataframe, if we are running the iteration process, we incorporate that info into the file name
    if iteration == False:
        mean_GnomAD_matrix.to_csv(f"{cwd_path}Motifs/{motif_name}/MeanGnomADMatrix_{motif_name}_positives.tsv", sep="\t")

        # plotting
        sns.set(font_scale=1.2)
        plt.figure(figsize=(15, 5))
        cmap = sns.color_palette("rocket", as_cmap=True)
        ax = sns.heatmap(mean_GnomAD_matrix, vmin=-6, vmax=-2, cmap=cmap, annot=True, linecolor="white", linewidths=.5, cbar_kws={'label':'-log10(Fo)'})
        plt.yticks(rotation=0)
        plt.ylabel("")
        plt.savefig(f"{cwd_path}Motifs/{motif_name}/MeanGnomADMatrix_HeatMap_{motif_name}_positives.png")
        if plot == True:
            plt.show()
        else:
            plt.close()
        
    else:
        mean_GnomAD_matrix.to_csv(f"{cwd_path}Motifs/{motif_name}/MeanGnomADMatrix_{motif_name}_positives_{iteration}_{number}.tsv", sep="\t")

        # plotting
        sns.set(font_scale=1.2)
        plt.figure(figsize=(15, 5))
        cmap = sns.color_palette("rocket", as_cmap=True)
        ax = sns.heatmap(mean_GnomAD_matrix, vmin=-6, vmax=-2, cmap=cmap, annot=True, linecolor="white", linewidths=.5, cbar_kws={'label':'-log10(Fo)'})
        plt.yticks(rotation=0)
        plt.ylabel("")
        plt.savefig(f"{cwd_path}Motifs/{motif_name}/MeanGnomADMatrix_HeatMap_{motif_name}_positives_{iteration}_{number}.png")
        if plot == True:
            plt.show()
        else:
            plt.close()

    return mean_GnomAD_matrix

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
    parser.add_argument("-gad_miss", "--gnomad_missense", action="store", help = "csv file with all gnomad missense variants", default="./GnomAD/GnomAD_missense.csv")
    parser.add_argument("-af_thr", "--af_threshold", action="store", help = "csv file with allele frequency for each UniProt, SwissProt entry gene", default="./GnomAD/uniprot_genes_freq_cutoffs.csv")
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

    from GnomADMotif import GnomADMotif
    GnomADMotif(screening_proteoma, args.Motif_name, args.gnomad_missense, args.af_threshold)

    GnomADMatrix = get_GnomADMatrix(screening_proteoma, args.Motif_label, args.Motif_name, args.cwd)
