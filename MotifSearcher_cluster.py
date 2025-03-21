#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
from Bio import SeqIO
import re
import argparse
from tqdm import tqdm
import pandas as pd
import os

class motif_object(object):
    """
    Object with all the match information for a specific motif in a given sequence of the human proteome.
    Important features: UniProtID of the protein, specific sequence of the motif in the match, start and end
    of the motif in one-based numbering, length of the protein, list with all strings used for gene definition.
    There are attributes born with default values that we are going to modify with methods in the future
    if certain criteria are met.
    """
    def __init__(self, UniProtID, hit, xml_record, UniProtIDs_omim_pheno):
        self.UniProtID = UniProtID
        self.seq = hit.group()
        self.start = hit.start() + 1 # numeración base 1
        self.end = hit.end() + 1 # numeración base 1, un aa más que el final del motif (listo para range)
        self.length = len(xml_record.seq)
        if UniProtID in UniProtIDs_omim_pheno:
            self.OMIM = True
        else:
            self.OMIM = False

        self.gene = []
        if "gene_name_primary" in xml_record.annotations.keys():
            self.gene.append(xml_record.annotations["gene_name_primary"])
        if "gene_name_synonym" in xml_record.annotations.keys():
            for gene in xml_record.annotations["gene_name_synonym"]:
                self.gene.append(gene)

        self.ELM_Positive = None
        self.ELM_Negative = None
        self.PDBs = []
        self.MotSASi_Positive = None
        self.MotSASi_Negative = None
        self.conflict = False
        self.passed_matrix = None
        self.passed_conservation = None
        self.passed_secondary_structure = None
        self.passed_relative_position = None
        self.passed_GOTerms = None
        self.passed_exposure = None
        self.Variants = []
        self.number_MotSASi_variants = len([variant for variant in self.Variants if variant.verdict in ["Benign", "Pathogenic"]])

    def call_true_positive(self):
        self.ELM_Positive = True

    def call_false_positive(self):
        self.ELM_Negative = True

    def call_MotSASi_positive(self):
        self.MotSASi_Positive = True

    def call_MotSASi_negative(self):
        self.MotSASi_Negative = True

    def add_pdbs(self, pdb):
        self.PDBs.append(pdb)

    def add_variants(self, variant):
        self.Variants.append(variant)
        self.number_MotSASi_variants = len([variant for variant in self.Variants if variant.verdict in ["Benign", "Pathogenic"]])

    def conservation_label(self, label):
        self.passed_conservation = label

    def secondary_structure_label(self, label):
        self.passed_secondary_structure = label

    def relative_position_label(self, label):
        self.passed_relative_position = label

    def GOTerms_label(self, label):
        self.passed_GOTerms = label

    def exposure_label(self, label):
        self.passed_exposure = label
    
    def conflict_label(self, label):
        self.conflict = label

    def matrix_label(self, label):
        self.passed_matrix = label

    def revert_MotSASi_verdict(self):
        self.MotSASi_Positive = False
        self.MotSASi_Negative = True

    def __eq__(self, other) :
        return self.__dict__ == other.__dict__

def get_acc_sp(record):
    """
    Esta funcion permite que al indexar SwissProt las keys del diccionario sean los UniProtID
    """
    parts = record.split("|")
    assert len(parts) == 3 and parts[0] == "sp"
    return parts[1]

def HumanProteomeSearch(motif_name, motif_re, sp_h_fasta_path, sp_h_xml_path, sp_h_list_path, uniprot_omim_phen_path):
    """
    Searches for motif regular expression in the human proteome. Outputs the hits in
    a list of objects with information regarding the protein that contains
    the motif, its location, etc. Important information is printed.
    In: motif name, motif regular expression, UniProt Swiss-Prot xml database and list with UniProtIDs on it.
    Out: list with motif objects.
    """

    print("Looking for coincidences of the motif in the human proteome")

    # parse UniProt Swiss-Prot fasta and xml database and generate the list with UniProtIDs contained on it
    fastas = SeqIO.index_db(f"{sp_h_fasta_path}.idx", f"{sp_h_fasta_path}", "fasta", key_function=get_acc_sp)
    xmls = SeqIO.index_db(f"{sp_h_xml_path}.idx", f"{sp_h_xml_path}", "uniprot-xml")
    proteoma_sp = pd.read_csv(f"{sp_h_list_path}", header=None)
    lista_proteoma_sp = list(proteoma_sp[0])
    uniprot_omim_phen = pd.read_csv(uniprot_omim_phen_path)
    UniProtIDs_omim_pheno = uniprot_omim_phen.UniProtID.to_list()

    # list with all proteome matches (as objects)
    screening_proteoma = []

    # iterate through the proteome and load the fasta entry (faster)
    for UniProtID in tqdm(lista_proteoma_sp):
        fasta_record = fastas[UniProtID]
        # search of the regular expression on the entry sequence
        for hit in re.finditer(motif_re, str(fasta_record.seq)):
            # for each coincidence we generate a match object and add it to the previous list, 
            # only if there is at least a gene associated
            # we use the xml entry now because we are interested in gene primary names and synonyms
            xml_record = xmls[UniProtID]
            match = motif_object(UniProtID, hit, xml_record, UniProtIDs_omim_pheno)
            if match.gene != []:
                screening_proteoma.append(match)

    print(f"{motif_name} motif found in {str(len(set([match.UniProtID for match in screening_proteoma])))} human proteins")
    print(f"Total hits in the human proteome: {str(len(screening_proteoma))}")

    return screening_proteoma

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Returns a list of match objects for the motif on a sequences database (default SwissProt)")
    parser.add_argument("Motif_name", action="store", help = "Name of the motif, e.g. DOC_MAPK_JIP1_4")
    parser.add_argument("Motif_re", action="store", help = "Motif regular expression, e.g. [RK]P[^P][^P]L.[LIVMF]")
    parser.add_argument("Motif_label", action="store", help = "Motif regular expression for plotting, e.g. RK.P.^P.^P.L.x.LIVMF")
    parser.add_argument("-fasta", "--sp_h_fasta", action="store", help = "file with fasta entries for SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.fasta")
    parser.add_argument("-xml", "--sp_h_xml", action="store", help = "file with xml entries for SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.xml")
    parser.add_argument("-list", "--sp_h_list", action="store", help = "file with all UniProtIDs in SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.list")
    parser.add_argument("-tmp", "--tmp_path", action="store", help = "temporary folder", default="/grupos/Marce/estructural/motsasi/tmp/")
    parser.add_argument("-cwd", "--cwd", action="store", help = "current working directory path", default="/grupos/Marce/estructural/motsasi/")

    args = parser.parse_args()
    
    if not os.path.exists(f"{args.cwd}Motifs/{args.Motif_name}"):
        os.system(f"mkdir {args.cwd}Motifs/{args.Motif_name}")

    if not os.path.exists(f"{args.tmp_path}"):
        os.system(f"mkdir {args.tmp_path}")

    screening_proteoma = HumanProteomeSearch(args.Motif_name, args.Motif_re, args.sp_h_fasta, args.sp_h_xml, args.sp_h_list)
