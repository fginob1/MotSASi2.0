#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import os
import urllib
import matplotlib as mpl
from tqdm import tqdm
import argparse

mpl.rcParams.update(mpl.rcParamsDefault)

def get_acc_sp(record):
    """
    This function allows to index SwissProt database with UniProtIDs as list
    """
    parts = record.split("|")
    assert len(parts) == 3 and parts[0] == "sp"
    return parts[1]

def get_acc_tr(record):
    """
    This function allows to index Trembl database with UniProtIDs as list
    """
    parts = record.split("|")
    assert len(parts) == 3 and parts[0] == "tr"
    return parts[1]

def get_bitscores_list(match, orthologous_seq, M_length, cwd_path, tmp_path):
    """
    This function generates a MSA from the orthologous sequences and computes conservation bitscores, then
    outputs a list with those bitscores in a list with values of the corresponding positions
    """
    # write a fasta file from the SeqRecord objects
    SeqIO.write(orthologous_seq, f"{tmp_path}archivo_seqs_{match.UniProtID}.fasta", "fasta")
    # generate the MSA
    os.system(f"mafft --quiet {tmp_path}archivo_seqs_{match.UniProtID}.fasta > {tmp_path}mafft_{match.UniProtID}.fasta")
    # run the algorithm to generate conservation bitscores for the match motif residues, 
    # output is a tsv file that we parse as a pandas dataframe
    os.system(f"python3 {cwd_path}Motsasi_jsdiv.py {match.UniProtID} {tmp_path}")
    js = pd.read_csv(f"{tmp_path}Conservacion_jsdiv_{match.UniProtID}.tsv", header=0, sep="\t", index_col="Residue_number")
    # generate a np.nan list with equal length from that of the motif
    lista_aux = [np.nan]*M_length
    # iterate through motif residues
    for n, aa in enumerate(match.seq):
        # keep the bitscore of the residue
        bitscore = js.loc[match.start + n].Bitscores
        # check if the residue is in a poorly represented position of the MSA (lots of gaps, returns -1000)
        if bitscore != -1000:
            # replace the nan value in the motif position index by the bitscore
            lista_aux[n] = bitscore
        else:
            continue
    
    return lista_aux

def get_orthologous_seq_list(match, UniRef, lista_sp, sp, trembl):
    """
    This function gets all the orthologous sequences from UniRef90 as a list of SeqRecord objects,
    if there is mo UniRef90 Cluster, return False
    """
    # filter UniRef90 database by UniProtID
    UniRef_UniProtID = UniRef[UniRef["UniProtID"] == match.UniProtID]
    # if find something go on, else print a message and continue
    if not UniRef_UniProtID.empty:
        # take the cluster participants, split them, and keep those that does not start start with UPI
        # (they belong to UniParc)
        ortologos_ids = UniRef_UniProtID["UniRef50_Cluster"].iloc[0].split(",")
        ortologos_ids.insert(0, match.UniProtID)
        # list to keep the SeqRecord objects with the working sequences
        ortologos_seq = []
        # list to check that we only keep one sequence by species
        ortologos_organisms = []
        # iterate through the UniProtIDs of the cluster
        for ortologo_id in ortologos_ids:
            # look for the UniProtID fasta format entry, if it is not present in SwissProt, go to Trembl
            if ortologo_id in lista_sp:
                seq_record = sp[ortologo_id]
            else:
                # must use try because many times the cluster member has been eliminated from Trembl
                try:
                    seq_record = trembl[ortologo_id]
                except:
                    continue
            # want to know species, we split the description and look for it
            desc_split = seq_record.description.split("=")
            organism = [desc_split[i+1][:-3] for i, n in enumerate(desc_split) if n.endswith("OS")][0]
            # check if had already incorporated a sequence from that species, check for length threshold too
            #if organism not in ortologos_organisms and len(seq_record.seq) > match.length*0.1:
            if organism not in ortologos_organisms:
                # append the SeqRecord object to the ortologos_seq list, id equal to UniProtID
                ortologos_seq.append(SeqRecord(Seq(str(seq_record.seq)), id=ortologo_id, description=""))
                # if append the sequence, append the species to the other list too
                ortologos_organisms.append(organism)

        return ortologos_seq

    else:
        return False

def ConservationScore(matchs_list, motif_label, motif_name, uniref50_db_path, sp_list_path, sp_fasta_path, trembl_fasta_path, cwd_path, tmp_path, plot = True, screen = False, save_tsv = True):
    """
    This function takes positive control motifs and generates a MSA doing use of their UniRef90 entries,
    with the goal of getting a conservation score by means of using the Jensen-Shannon Divergence algorithm
    """

    # length of the motif
    M_length = len(motif_label.split("."))

    # load the UniRef90 database as a pandas dataframe and clean it to have accessible UniProtIDs
    UniRef = pd.read_csv(uniref50_db_path)
    #UniRef["Cluster ID"].replace("UniRef90_","", regex=True, inplace=True)

    # load the UniProtIDs from SwissProt as a list
    raw_sp = pd.read_csv(sp_list_path, header=None)
    lista_sp = set(raw_sp[0])

    # parse the SwissProt database in fasta format
    sp = SeqIO.index_db(f"{sp_fasta_path}.idx", sp_fasta_path, "fasta", key_function=get_acc_sp)

    # parse the Trembl database in fasta format
    trembl = SeqIO.index_db(f"{trembl_fasta_path}.idx", trembl_fasta_path, "fasta", key_function=get_acc_tr)

    # dictionary to keep conservation score results, keys are UniProtIDs, values are lists with scores for every position
    score = {}
    
    if len(matchs_list) == 1:
        # simple way to not modify the code that much, it is just a match
        for match in matchs_list:
            # get a orthologous list of sequences from different species
            orthologous_seq = get_orthologous_seq_list(match, UniRef, lista_sp, sp, trembl)
            if orthologous_seq == False:
                # continue if there is no UniRef90 Cluster
                print(f"Motif found in {match.UniProtID} has no UniRef90 cluster, it will not be taken into consideration")
                continue
            # check if there are at least two sequences
            if len(orthologous_seq) > 2:
                # get a list with length equal to the motif length with conservation bitscore to append it
                # to the dictionary with UniProtID of the march as key
                lista_aux = get_bitscores_list(match, orthologous_seq, M_length, cwd_path, tmp_path)
                # keep the info in the score dictionary if all positions contribute with bitscores
                if np.nan not in lista_aux:    
                    score[match.UniProtID] = lista_aux
            else:
                print(f"Motif found in {match.UniProtID} presents 2 or less orthologs, it will not be taken into consideration")

    else:
        # iterate through positive control motifs
        for match in tqdm(matchs_list):
            orthologous_seq = get_orthologous_seq_list(match, UniRef, lista_sp, sp, trembl)
            if orthologous_seq == False:
                # continue if there is no UniRef90 Cluster
                print(f"Motif found in {match.UniProtID} has no UniRef90 cluster, it will not be taken into consideration")
                continue
            # check if there are at least two sequences
            if len(orthologous_seq) > 2:
                # get a list with length equal to the motif length with conservation bitscore to append it
                # to the dictionary with UniProtID of the march as key
                lista_aux = get_bitscores_list(match, orthologous_seq, M_length, cwd_path, tmp_path)
                # keep the info in the score dictionary if all positions contribute with bitscores
                if np.nan not in lista_aux:    
                    score[match.UniProtID] = lista_aux
            else:
                print(f"Motif found in {match.UniProtID} presents 2 or less orthologs, it will not be taken into consideration")

    # generate a pandas dataframe from the score dictionary
    df_score = pd.DataFrame(score)
    # rows labels, create a column and set them as index
    df_score["Motif"] = motif_label.split(".")
    df_score.set_index("Motif", inplace=True, drop=True)

    if save_tsv == True:
        df_score.to_csv(f"{cwd_path}Motifs/{motif_name}/ConservationScore_{motif_name}_positives.tsv", sep="\t")

    if plot == True:
        # must generate an index with no repetitions in Motifs, otherwise plot construction will fail
        # solution by adding the index position
        df_score["Motif"] = [res +"_"+ str(n) for n, res in enumerate(motif_label.split("."))]
        df_score.set_index("Motif", inplace=True, drop=True)
        # plotting
        plt.figure(figsize=(8, 4), facecolor="w")
        ax = sns.boxplot(data=df_score.T, color="white", width=0.2, orient="h")
        plt.xlabel("Conservation Score")
        plt.xlim(0.65, 0.95)
        ax.set_yticklabels(motif_label.split("."))
        plt.savefig(f"{cwd_path}Motifs/{motif_name}/ConservationScore_{motif_name}_positives.png", bbox_inches = "tight")
        if screen == True:
            plt.show()
        else:
            plt.close()

        # recover the Motif index
        df_score["Motif"] = motif_label.split(".")
        df_score.set_index("Motif", inplace=True, drop=True)

    return df_score

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
    parser.add_argument("-uniref90", "--uniref90_clusters", action="store", help = "file with UniRef90 clusters for each UniProtID", default="./UniProt/UniRef-HomoSapiens-90.tsv")
    parser.add_argument("-sp_list", "--sp_list", action="store", help = "file with all UniProtIDs in SwissProt", default="./UniProt/uniprot_sprot.list")
    parser.add_argument("-sp_fasta", "--sp_fasta", action="store", help = "file with all fasta entries in SwissProt", default="./UniProt/uniprot_sprot.fasta.bgz")
    parser.add_argument("-trembl", "--trembl_fasta", action="store", help = "file with all fasta entries in TrEMBL", default="./UniProt/uniprot_trembl.fasta.bgz")
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
    print("Scanning Motif Conservation")
    ConservationScore(true_positives, args.Motif_label, args.Motif_name, args.uniref90_clusters, args.sp_list, args.sp_fasta, args.trembl_fasta, args.cwd, args.tmp_path)
