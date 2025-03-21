#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
from itertools import repeat
import argparse
import os

class variante_clinvar(object):
    """
    The ClinVar variant object carry with the most important information about variants:
    UniProtID, reference and alternative residue, position, clinical significance,
    review status, accession, position in the motif.
    The verdict attribute is only with purpose of classify variants as pathogenic or benign and source_db
    has the only purpose of knowing the source
    """
    def __init__(self, UniProtID, variant, verdict, j, motif_seq, ELM_Positive):
        self.UniProtID = UniProtID
        self.ref = variant.Referencia
        self.pos = variant.Posicion
        self.alt = variant.Alternativo
        self.clinsig = variant.Clinsig
        self.verdict = verdict
        self.accession = variant.Accession
        self.source_db = "ClinVar"
        self.index = j
        self.ddGs_binding = []
        self.ddGs_stability = []
        self.motif_seq = motif_seq
        self.motif_ELM_Positive = ELM_Positive
        self.status = variant.ReviewStatus
        self.submissions = variant.Number_of_submissions
        confidence_score = self.submissions * 0.1
        if self.status in ["no assertion provided", "no assertion criteria provided", "no assertion for the individual variant"]:
            confidence_score += 0
        elif self.status in ["criteria provided, single submitter", "criteria provided, conflicting interpretations"]:
            confidence_score += 1
        elif self.status in ["criteria provided, multiple submitters, no conflicts"]:
            confidence_score += 2
        elif self.status in ["reviewed by expert panel"]:
            confidence_score += 3
        elif self.status in ["practice guideline"]:
            confidence_score += 4
        if self.verdict == "Pathogenic":
            self.confidence = confidence_score
        elif self.verdict == "Benign":
            self.confidence = -confidence_score

    def add_ddGb(self, ddGb):
        self.ddGs_binding.append(ddGb)

    def add_ddGs(self, ddGs):
        self.ddGs_stability.append(ddGs)


def clinvar_search(UniProtID, screening_proteoma, var_ClinVar):    
    """
    This function looks for ClinVar variants ocurring within the motif of a particular protein.
    Once found, the variant changs its shape to an object (ClinVar variant object) and we append it to the motif object.
    In: UniProtID, list of motif objects, file of ClinVar missense variants.
    Out: list with ClinVar variant objects.
    """

    # look for the list of genes associated to that UniProtID
    for match in screening_proteoma:
        if match.UniProtID == UniProtID:
            gene = match.gene
            break

    # matchs in screening_proteoma that correspond to that UniProtID
    matchs_UniProtID = [match for match in screening_proteoma if match.UniProtID == UniProtID]

    # filter the ClinVar file by variants located on that gene (remember it is a list with different ways of naming the gene)
    var_ClinVar_gene = var_ClinVar[var_ClinVar["Gene"].isin(gene)]
    if not var_ClinVar_gene.empty:
        # list where we are going to append the ClinVar variant objects
        variantes_clinvar = []
        # iterate through matches, filter the variants by position in the protein and, if there is a coincidence
        # of position and reference residue, we generate the ClinVar variant object and both append it to the list
        # and add it to the variants list of the match
        for match in matchs_UniProtID:
            var_ClinVar_gene_pos = var_ClinVar_gene[(var_ClinVar_gene.Posicion >= match.start) & (var_ClinVar_gene.Posicion < match.end)]
            if not var_ClinVar_gene_pos.empty:
                for i, variant in var_ClinVar_gene_pos.iterrows():
                    for j, aa in enumerate(match.seq):
                        if variant.Referencia == aa and variant.Posicion == match.start + j:
                            if variant.Clinsig in ["Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"]:
                                variantes_clinvar.append(variante_clinvar(UniProtID, variant, "Pathogenic", j, match.seq, match.ELM_Positive))
                                match.add_variants(variante_clinvar(UniProtID, variant, "Pathogenic", j, match.seq, match.ELM_Positive))
                            elif variant.Clinsig in ["Benign", "Benign/Likely benign", "Likely benign"]:
                                variantes_clinvar.append(variante_clinvar(UniProtID, variant, "Benign", j, match.seq, match.ELM_Positive))
                                match.add_variants(variante_clinvar(UniProtID, variant, "Benign", j, match.seq, match.ELM_Positive))
                            else:
                                match.add_variants(variante_clinvar(UniProtID, variant, "Uncertain", j, match.seq, match.ELM_Positive))

        return [variantes_clinvar, matchs_UniProtID]

    else:
        return [[], matchs_UniProtID]


def ClinVarMotif(screening_proteoma, motif_name, clinvar_missense_bp_path):
    """
    List of ClinVar variant objects that occur within the motif
    """

    print("Inspecting ClinVar...")

    # parsing the clinvar missense database previously generated
    var_ClinVar = pd.read_csv(clinvar_missense_bp_path)

    # list with unique UniProtIDs in the motif objects list
    UniProtIDs = list(set([match.UniProtID for match in screening_proteoma]))

    # list with all ClinVar variants on studied potential motifs
    ClinVar_variant_objects = []
        
    with Pool() as pool:
        output_ClinVar = pool.starmap(clinvar_search, zip(UniProtIDs, repeat(screening_proteoma), repeat(var_ClinVar)))
    
    ClinVar_variant_objects = [variant for aux_list_1 in output_ClinVar for variant in aux_list_1[0]]

    screening_proteoma = [match for aux_list_1 in output_ClinVar for match in aux_list_1[1]]

    if len(ClinVar_variant_objects) > 0:
        bp__clinvar_variants = len([variant for match in screening_proteoma for variant in match.Variants if variant.source_db == "ClinVar" and variant.verdict in ["Benign", "Pathogenic"]])
        total_clinvar_variants = len([variant for match in screening_proteoma for variant in match.Variants if variant.source_db == "ClinVar"])
        print(f"Total number of Benign/Pathogenic ClinVar variants affecting {motif_name} motif: {bp__clinvar_variants}")
        print(f"Total number of ClinVar variants affecting {motif_name} motif: {total_clinvar_variants}")
    else:
        print(f"There are no Benign/Pathogenic ClinVar variants reported within these motifs...")

    return ClinVar_variant_objects, screening_proteoma

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Returns a list of match objects for the motif on a sequences database (default SwissProt)")
    parser.add_argument("Motif_name", action="store", help = "Name of the motif, e.g. DOC_MAPK_JIP1_4")
    parser.add_argument("Motif_re", action="store", help = "Motif regular expression, e.g. [RK]P[^P][^P]L.[LIVMF]")
    parser.add_argument("Motif_label", action="store", help = "Motif regular expression for plotting, e.g. RK.P.^P.^P.L.x.LIVMF")
    parser.add_argument("-fasta", "--sp_h_fasta", action="store", help = "file with fasta entries for SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.fasta")
    parser.add_argument("-xml", "--sp_h_xml", action="store", help = "file with xml entries for SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.xml")
    parser.add_argument("-list", "--sp_h_list", action="store", help = "file with all UniProtIDs in SwissProt", default="/grupos/Marce/estructural/databases/UniProt/uniprot_sprot_h.list")
    parser.add_argument("-cv_miss", "--clinvar_missense", action="store", help = "csv file with all clinvar missense variants", default="/grupos/Marce/estructural/databases/ClinVar/ClinVar_missense_all_filtered.csv")
    parser.add_argument("-tmp", "--tmp_path", action="store", help = "temporary folder", default="/grupos/Marce/estructural/motsasi/tmp/")
    parser.add_argument("-cwd", "--cwd", action="store", help = "current working directory path", default="/grupos/Marce/estructural/motsasi/")

    args = parser.parse_args()
    
    if not os.path.exists(f"{args.cwd}Motifs/{args.Motif_name}"):
        os.system(f"mkdir {args.cwd}Motifs/{args.Motif_name}")

    if not os.path.exists(f"{args.tmp_path}"):
        os.system(f"mkdir {args.tmp_path}")

    from MotifSearcher_cluster import HumanProteomeSearch
    screening_proteoma = HumanProteomeSearch(args.Motif_name, args.Motif_re, args.sp_h_fasta, args.sp_h_xml, args.sp_h_list)

    ClinVarMotif(screening_proteoma, args.Motif_name, args.clinvar_missense)
