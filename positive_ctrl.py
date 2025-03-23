#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
from Bio.PDB import PDBParser
import pandas as pd
import gzip
import os
import warnings
import argparse
warnings.filterwarnings("ignore")

def positive_ctrl(screening_proteoma, motif_name, pdb_db_path, elm_instances_path, overlap = 10):
    """
    Looks for the positive control group of the Motif.
    These are the motifs that have been experimentally tested.
    In: list with motif objects, motif name (as ELMIdentifier), path to pdb crystal structures, path to elm_instances database
    Out: nothing, and pdb info are labeled or added to the motif objects while parsing them
    """
    # parsing elm_instances database and filter by motif name and organism
    elm = pd.read_csv(elm_instances_path, index_col=None, skiprows=5, sep="\t")
    elm_filtered = elm[(elm["ELMIdentifier"] == motif_name) & (elm["Organism"] == "Homo sapiens")]
    # keep UniProtIDs (Primary_Acc) of matches, useful for work, that we can find in our motif objects list, without modifying them
    UniProtIDs = [x for x in elm_filtered["Primary_Acc"].tolist() if "-" not in x or x.endswith("-1")]
    UniProtIDs_filtered = [x for x in UniProtIDs if x.split("-")[0] in set([match.UniProtID for match in screening_proteoma])]
    # filter by those UniProtIDs
    elm_filtered = elm_filtered[elm_filtered["Primary_Acc"].isin(UniProtIDs_filtered)]
    # clean UniProtIDs
    elm_filtered["Primary_Acc"] = [x.split("-")[0] for x in elm_filtered["Primary_Acc"].tolist()]

    # the database contains true positives (most) as false positives too (minority), we keep them by separate
    true_positives = elm_filtered[elm_filtered.InstanceLogic == "true positive"]

    false_positives = elm_filtered[elm_filtered.InstanceLogic == "false positive"]

    # iterate through our motif objects list
    for match in screening_proteoma:
        # iterate through true positives in order to try to label the motif object as a positive control
        # done that we have seen disagreements in numbering, we include the overlap parameter that allows us
        # to keep the positive control if the difference in position numbering is small
        for i, row in true_positives.iterrows():
            if match.UniProtID == row.Primary_Acc and row.Start-overlap <= match.start <= row.Start+overlap:
                match.call_true_positive()
        # same for false positives
        for i, row in false_positives.iterrows():
            if match.UniProtID == row.Primary_Acc and row.Start-overlap <= match.start <= row.Start+overlap:
                match.call_false_positive()

    # filter for matches with pdb info and iterate through them checking if those crystal structures correspond
    # to x-ray diffraction experiments (recommended for FoldX calculations)
    elm_pdb = elm_filtered[elm_filtered.PDB.notnull()]
    for i, row in elm_pdb.iterrows():
        pdbs = row.PDB.split()
        for pdb in pdbs:
            parser = PDBParser(PERMISSIVE=1)
            structure_id = pdb
            filename = f"{pdb_db_path}pdb{pdb.lower()}.ent.gz"
            structure = parser.get_structure(structure_id, gzip.open(filename, "rt"))
            if structure.header["structure_method"] == "x-ray diffraction":
                for match in [match for match in screening_proteoma if match.ELM_Positive == True]:
                    if match.UniProtID == row.Primary_Acc and row.Start-overlap <= match.start <= row.Start+overlap and pdb not in match.PDBs:
                        # we append the pdb identifier to a list attribute of the specific match
                        match.add_pdbs(pdb)

    print(f"ELM True positives: {str(len([match for match in screening_proteoma if match.ELM_Positive == True]))}")
    print(f"ELM True negatives: {str(len([match for match in screening_proteoma if match.ELM_Negative == True]))}")
    print(f"Matchs containing PDBs associated: {str(len([match for match in screening_proteoma if len(match.PDBs) > 0]))}")

    return

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
    parser.add_argument("-tmp", "--tmp_path", action="store", help = "temporary folder", default="/grupos/Marce/estructural/motsasi/tmp/")
    parser.add_argument("-cwd", "--cwd", action="store", help = "current working directory path", default="/grupos/Marce/estructural/motsasi/")

    args = parser.parse_args()

    if not os.path.exists(f"{args.cwd}Motifs/{args.Motif_name}"):
        os.system(f"mkdir {args.cwd}Motifs/{args.Motif_name}")

    if not os.path.exists(f"{args.tmp_path}"):
        os.system(f"mkdir {args.tmp_path}")

    from MotifSearcher import HumanProteomeSearch
    screening_proteoma = HumanProteomeSearch(args.Motif_name, args.Motif_re, args.sp_h_fasta, args.sp_h_xml, args.sp_h_list)

    positive_ctrl(screening_proteoma, args.Motif_name, args.pdb_structures, args.elm_instances)
