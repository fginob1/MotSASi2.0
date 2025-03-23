#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
import pandas as pd
import math
from tqdm import tqdm
import os
from multiprocessing import Pool
from itertools import repeat
import argparse

class variante_gnomad(object):
    """
    The GnomAD variant object carry with the most important information about variants:
    UniProtID, reference and alternative residue, position, allele frequency, position in the motif.
    The verdict attribute is only with purpose of classify variants as pathogenic or benign and source_db
    has the only purpose of knowing the source
    """
    def __init__(self, UniProtID, variant, verdict, j, switch, min_AF, motif_seq, ELM_Positive):
        self.UniProtID = UniProtID
        self.ref = variant.Referencia
        self.pos = variant.Posicion
        self.alt = variant.Alternativo
        self.freq = variant.AF
        self.max_freq = variant.AF_popmax
        #self.minus_log_max_freq = -math.log10(variant.AF_popmax)
        #self.confidence = (1/3)*(math.sqrt(3))**(-math.log10(variant.AF_popmax)) - 5
        self.minus_log_max_freq = -math.log10(variant.AF)
        self.confidence = (1/3)*(math.sqrt(3))**(-math.log10(variant.AF)) - 5
        self.verdict = verdict
        self.source_db = "GnomAD"
        self.index = j
        if switch == True:
            self.freq = min_AF
        self.ddGs_binding = []
        self.ddGs_stability = []
        self.motif_seq = motif_seq
        self.motif_ELM_Positive = ELM_Positive

    def add_ddGb(self, ddGb):
        self.ddGs_binding.append(ddGb)

    def add_ddGs(self, ddGs):
        self.ddGs_stability.append(ddGs)


def gnomad_search(UniProtID, screening_proteoma, var_GnomAD, sp_h_genes_af_cut_off_path):
    """
    This function looks for GnomAD variants ocurring within the motif of a particular protein.
    Once found, the variant changs its shape to an object (GnomAD variant object) and we append it to the motif object.
    In: UniProtID, list of motif objects, file of GnomAD missense variants.
    Out: list with ClinVar variant objects.
    """
    # get the allele frequency threshold for the corresponding UniProtID
    af_cut_off_df = pd.read_csv(sp_h_genes_af_cut_off_path, index_col=["UniProtID"])
    af_thr = af_cut_off_df.loc[UniProtID, "Freq_cut_off"]

#     # look for the list of genes associated to that UniProtID
#     for match in screening_proteoma:
#         if match.UniProtID == UniProtID:
#             gene = match.gene
#             break
#     # filter the GnomAD file by variants located on that gene (remember it is a list with different ways of naming the gene)
#     var_GnomAD_gene = var_GnomAD[var_GnomAD["Gene"].isin(gene)]

    # filter the GnomAD file by variants located on that gene (remember it is a list with different ways of naming the gene)
    var_GnomAD_gene = var_GnomAD[var_GnomAD["UniProtID"] == UniProtID]

    # minimum allele frequency detected on the database for that protein missense variants
    min_AF = var_GnomAD_gene.AF_popmax[var_GnomAD_gene.AF_popmax > 0].min()
    #min_AF = var_GnomAD_gene.AF[var_GnomAD_gene.AF > 0].min()
    
    # matchs in screening_proteoma that correspond to that UniProtID
    matchs_UniProtID = [match for match in screening_proteoma if match.UniProtID == UniProtID]

    if not var_GnomAD_gene.empty:
        # list where we are going to append the GnomAD variant objects (always benign)
        benign_gnomad_variants = []
        # iterate through matches, filter the variants by position in the protein and, if there is a coincidence
        # of position and reference residue, we generate the GnomAD variant object and both append it to the list
        # and add it to the variants list of the match, only if the threshold is met
        # if the AF is equal to 0, we keep the minimum allele frequency
        for match in matchs_UniProtID:
            var_GnomAD_gene_pos = var_GnomAD_gene[(var_GnomAD_gene.Posicion >= match.start) & (var_GnomAD_gene.Posicion < match.end)]
            if not var_GnomAD_gene_pos.empty:
                for i, variant in var_GnomAD_gene_pos.iterrows():
                    for j, aa in enumerate(match.seq):
                        if variant.Referencia == aa and variant.Posicion == match.start + j:
                            #if variant.AF_popmax != 0:
                            if variant.AF != 0:
                                #if -math.log10(variant.AF_popmax) <= af_thr:
                                if -math.log10(variant.AF) <= af_thr:
                                    #print(match.__dict__)
                                    benign_gnomad_variants.append(variante_gnomad(UniProtID, variant, "Benign", j, False, min_AF, match.seq, match.ELM_Positive))
                                    match.add_variants(variante_gnomad(UniProtID, variant, "Benign", j, False, min_AF, match.seq, match.ELM_Positive))
                                    #print(variante_gnomad(UniProtID, variant, "Benign", j, False, min_AF, match.seq).__dict__)
                                else:
                                    match.add_variants(variante_gnomad(UniProtID, variant, "Uncertain", j, False, min_AF, match.seq, match.ELM_Positive))
                            else:
                                match.add_variants(variante_gnomad(UniProtID, variant, "Uncertain", j, True, min_AF, match.seq, match.ELM_Positive))

        return [benign_gnomad_variants, matchs_UniProtID]

    else:
        return [[], matchs_UniProtID]

def GnomADMotif(screening_proteoma, motif_name, gnomad_missense_path, sp_h_genes_af_cut_off_path):
    """
    Iterates over chromosomes and creates the final list of GnomAD variants
    with -log10(AF)<4 that occurr within the motif
    """

    print("Inspecting GnomAD...")

    # list with unique UniProtIDs in the motif objects list
    UniProtIDs = list(set([match.UniProtID for match in screening_proteoma]))

    # parsing the gnomad missense database previously generated
    var_GnomAD = pd.read_csv(gnomad_missense_path)

    # filtering the gnomad missense database by the UniProtIDs of interest
    var_GnomAD = var_GnomAD[var_GnomAD.UniProtID.isin(UniProtIDs)]

    # list with all ClinVar variants on studied potential motifs
    #GnomAD_variant_objects = []

    # iterate through UniProtIDs to look for GnomAD variants and append them to the list as GnomAD variant objects
#     for UniProtID in tqdm(UniProtIDs):
#         GnomAD_variant_objects += gnomad_search(UniProtID, screening_proteoma, var_GnomAD, sp_h_genes_af_cut_off_path)

    with Pool() as pool:
        output_GnomAD = pool.starmap(gnomad_search, zip(UniProtIDs, repeat(screening_proteoma), repeat(var_GnomAD), repeat(sp_h_genes_af_cut_off_path)))

    GnomAD_variant_objects = [variant for aux_list_1 in output_GnomAD for variant in aux_list_1[0]]

    screening_proteoma = [match for aux_list_1 in output_GnomAD for match in aux_list_1[1]]

    if len(GnomAD_variant_objects) > 0:
        b_gnomad_variants = len([variant for match in screening_proteoma for variant in match.Variants if variant.source_db == "GnomAD" and variant.verdict in ["Benign"]])
        total_gnomad_variants = len([variant for match in screening_proteoma for variant in match.Variants if variant.source_db == "GnomAD"])
        print(f"Total number of benign GnomAD variants affecting {motif_name} motif: {b_gnomad_variants}")
        print(f"Total number of GnomAD variants affecting {motif_name} motif: {total_gnomad_variants}")
    else:
        print(f"There are no benign GnomAD variants reported within these motifs...")

    return GnomAD_variant_objects, screening_proteoma

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Returns a list of match objects for the motif on a sequences database (default SwissProt)")
    parser.add_argument("Motif_name", action="store", help = "Name of the motif, e.g. DOC_MAPK_JIP1_4")
    parser.add_argument("Motif_re", action="store", help = "Motif regular expression, e.g. [RK]P[^P][^P]L.[LIVMF]")
    parser.add_argument("Motif_label", action="store", help = "Motif regular expression for plotting, e.g. RK.P.^P.^P.L.x.LIVMF")
    parser.add_argument("-fasta", "--sp_h_fasta", action="store", help = "file with fasta entries for SwissProt", default="./UniProt/uniprot_sprot_h.fasta")
    parser.add_argument("-xml", "--sp_h_xml", action="store", help = "file with xml entries for SwissProt", default="./UniProt/uniprot_sprot_h.xml")
    parser.add_argument("-list", "--sp_h_list", action="store", help = "file with all UniProtIDs in SwissProt", default="./UniProt/uniprot_sprot_h.list")
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

    GnomADMotif(screening_proteoma, args.Motif_name, args.gnomad_missense, args.af_threshold)
