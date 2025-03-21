#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
import pandas as pd
from Bio import SeqIO
import collections
import requests
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
from tqdm import tqdm
import csv
import obonet

mpl.rcParams.update(mpl.rcParamsDefault)

def GO_terms_motif(UniProtID, xmls, GO_db, current_GO_Numbers):
    """
    This function returns all GO Terms associated to a given UniProtID (by UniProt)
    """

    # list to keep GO Terms
    GO_terms = []

    # parse the UniProt, SwissProt entry
    record = xmls[UniProtID]

    # iterate through database cross references and keep GO Numbers, then translate the Numbers to Terms
    # and append them to the list
    for entry in record.dbxrefs:
        if entry.startswith("GO:GO"):
            GO_Number = entry[3:]
            if GO_Number in current_GO_Numbers:
                if GO_db.nodes[GO_Number]["namespace"] == "cellular_component":
                    GO_terms.append(GO_db.nodes[GO_Number]["name"])

    return GO_terms


def Motif_GO(matchs_list, motif_name, sp_h_xml_path, GO_db_path, cwd_path, plot = True, screen = False, save_tsv = True):
    """
    This function identifies GO Terms associated to the positive control group (or a specific protein) 
    and allows to identify which of them are enriched (also by plotting)
    """

    # parse UniProt, SwissProt database in xml format
    xmls = SeqIO.index_db(f"{sp_h_xml_path}.idx", sp_h_xml_path, "uniprot-xml")
    
    # parse the GO database
    GO_db = obonet.read_obo(GO_db_path)
    current_GO_Numbers = set([term_id for term_id in GO_db.nodes])

    # list to keep GO Terms collected, count is a counter of proteins evaluated (note the set of UniProtIDs when len(matchs_list) > 1)
    GO_terms = []
    count = 0
    if len(matchs_list) == 1:
        # simple way to not modify the code that much, it is just an UniProtID
        for match in [match.UniProtID for match in matchs_list]:
            # call the function that outputs the GO Terms as a list and append them to the first list
            GO_terms += GO_terms_motif(match, xmls, GO_db, current_GO_Numbers)
            count += 1
    else:
        # iterate through the UniProtIDs of the positive control group
        for match in tqdm(set([match.UniProtID for match in matchs_list])):
            # call the function that outputs the GO Terms as a list and append them to the first list
            GO_terms += GO_terms_motif(match, xmls, GO_db, current_GO_Numbers)
            count += 1

    # generate a Counter object where GO Terms are counted
    GO_terms_counts = collections.Counter(GO_terms)
    # convert the Counter object into a dictionary where keys are ordered by value, which is the number of appearances
    GO_terms_counts = dict(sorted(GO_terms_counts.items(), key=lambda item: item[1], reverse=True))

    # generate a pandas dataframe from the dictionary, divide by the total number of proteins evaluated,
    # multiply by 100 to obtain percentages
    GO = pd.DataFrame(GO_terms_counts, index=["%"])
    GO = (GO/count)*100

    if save_tsv == True:
        GO.to_csv(f"{cwd_path}Motifs/{motif_name}/GOterms_{motif_name}_positives.tsv", sep="\t")

    if plot == True:
        # n allows to decide how many GO Terms we want when plotting
        n = 20
        if len(GO.columns) >= n:
            # plotting
            plt.figure()
            plt.title("GO terms - UniProt Keywords")
            plt.barh(GO.columns[:n], GO.iloc[0, :n])
            plt.yticks(range(n), GO.columns[:n], fontsize=10)
            plt.xlim((0,100))
            plt.xlabel("% of proteins of the positive control")
            plt.savefig(f"{cwd_path}Motifs/{motif_name}/GOterms_{motif_name}_positives.png", bbox_inches = "tight")
            if screen == True:
                plt.show()
            else:
                plt.close()

        else:
            # plotting
            plt.figure()
            plt.title("GO terms - UniProt Keywords")
            plt.barh(GO.columns, GO.iloc[0])
            plt.yticks(range(len(GO.columns)), GO.columns, fontsize=10)
            plt.xlim((0,100))
            plt.xlabel("% of proteins of the positive control")
            plt.savefig(f"{cwd_path}Motifs/{motif_name}/GOterms_{motif_name}_positives.png", bbox_inches = "tight")
            if screen == True:
                plt.show()
            else:
                plt.close()

    return GO

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
    parser.add_argument("-go_db", "--go_database", action="store", help = "file with codes for GO Terms", default="/grupos/Marce/estructural/databases/GO/go.obo")
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
    print("Scanning Motif GO Terms")
    Motif_GO(true_positives, args.Motif_name, args.sp_h_xml, args.go_database, args.cwd)
