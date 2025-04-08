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
from Bio.PDB import PDBParser
import numpy as np
import urllib
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import os
import matplotlib as mpl
from tqdm import tqdm
import ssbio.protein.sequence.properties.scratch
import collections
import requests
import json
import math
import gzip
import re
import csv
import argparse
from MotifSearcher import HumanProteomeSearch
from positive_ctrl import positive_ctrl
from ClinVarMotif import ClinVarMotif
from ClinVarMatrix import get_ClinVarMatrix
from GnomADMotif import GnomADMotif
from GnomADMatrix import get_GnomADMatrix
from FoldXMotif import foldx_process
from ConservationScore import ConservationScore
from SecStrMotif import scratch
from GOMotif import Motif_GO
from ExpBurMotif import sasa
from multiprocessing import Pool
from itertools import repeat

mpl.rcParams.update(mpl.rcParamsDefault)


class iteration_object_variant(object):
    """
    Objects used to keep variants info and source
    """
    def __init__(self, ref, index, alt, verdict):
        self.ref = ref
        self.index = index
        self.alt = alt
        self.verdict = verdict

    def __eq__(self, other) :
        return self.__dict__ == other.__dict__

def combine_dfs(dfs_list, aas):

    df_rows = []
    for n, (i, row) in enumerate(dfs_list[0].iterrows()):
        row_list = []
        for aa in dfs_list[0]:
            cell_list = [df.iloc[n][aa] for df in dfs_list]
            positive_values = sum([j for j in cell_list if j > 0])
            negative_values = sum([j for j in cell_list if j < 0])
            if abs(positive_values) > abs(negative_values):
                row_list.append(positive_values)
            else:
                row_list.append(negative_values)
        df_rows.append(row_list)

    df = pd.DataFrame(df_rows, columns=aas)

    return df

def P_B(motif_label, motif_name, motif_re, true_positives, positives, clinvar_db_path, gnomad_db_path, cwd_path, tmp_path, pdbs_reparados_path, pdb_db_path, iteration, number, screen = False):
    """
    This function takes both true positives from ELM and positives already analyzed by MotSASi and recovers
    all variants of interest from ClinVar and GnomAD on them
    """

    # aminoacids list in order
    aas = ["G", "A", "V", "L", "I", "M", "P", "W", "F", "Y", "S", "T", "C", "N", "Q", "D", "E", "K", "R", "H"]

    # list of indexes of the motif regular expression positions where we find a lot of degeneration
    flex_positions = [i for i, x in enumerate(motif_label.split(".")) if x.startswith("^") or x == "x"]

    # we eliminate duplicates from the sum of true positives and MotSASi positives
    true_and_new_positives = unique(true_positives + positives)
    print(f"True positives plus MotSASi positives: {len(true_and_new_positives)}")
#     for match in true_and_new_positives:
#         print(match.__dict__)

    # lists where putting variants from each search
    ClinVar_P = []
    ClinVar_B = []
    GnomAD_B = []
    FoldX_P = []
    FoldX_B = []

    # calculate the matrix from ClinVar variants
    ClinVarMatrix = get_ClinVarMatrix(true_and_new_positives, motif_label, motif_name, cwd_path, iteration = iteration, number = number)
    # we iterate over the ClinVar matrix and collect all pathogenic and benign variants into the lists
    for i in range(len(motif_label.split("."))):
        for alt in ClinVarMatrix:
            # Pathogenics
            if ClinVarMatrix.iloc[i][alt] > 0:
                # if we are in the more strict steps of the iteration process, we get variants only from
                # more conserved motif positions and only keep reference from regular expression residues
                if iteration in ["rigid", "foldx"] and i not in flex_positions:
                    motif_pos = motif_label.split(".")[i]
                    for aa in motif_pos:
                        ClinVar_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                # in the final iteration step we will keep variants no matter which residue was found in the
                # reference sequence
                elif iteration == "flexible":
                    for aa in aas:
                        ClinVar_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
            # Benigns
            elif ClinVarMatrix.iloc[i][alt] < 0:
                if iteration in ["rigid", "foldx"] and i not in flex_positions:
                    motif_pos = motif_label.split(".")[i]
                    for aa in motif_pos:
                        ClinVar_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                elif iteration == "flexible":
                    for aa in aas:
                        ClinVar_B.append(iteration_object_variant(aa, i, alt, "Benign"))

    # we simplify the numeric annotation getting just a three number code
    ClinVarMatrix.fillna(0, inplace=True)
    print(ClinVarMatrix)
    ClinVarMatrix.reset_index(inplace=True, drop=True)

    # calculate the matrix from GnomAD variants
    GnomADMatrix = get_GnomADMatrix(true_and_new_positives, motif_label, motif_name, cwd_path, iteration = iteration, number = number)
    # same for benign GnomAD variants
    for i in range(len(motif_label.split("."))):
        for alt in GnomADMatrix:
            # Benigns
            if not np.isnan(GnomADMatrix.iloc[i][alt]):
                if iteration in ["rigid", "foldx"] and i not in flex_positions:
                    motif_pos = motif_label.split(".")[i]
                    for aa in motif_pos:
                        GnomAD_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                elif iteration == "flexible":
                    for aa in aas:
                        GnomAD_B.append(iteration_object_variant(aa, i, alt, "Benign"))

    # same for GnomAD
    GnomADMatrix.fillna(0, inplace=True)
    print(GnomADMatrix)
    GnomADMatrix.reset_index(inplace=True, drop=True)

    # we sum both matrices
    FinalMatrix = combine_dfs([ClinVarMatrix, GnomADMatrix], aas)
    FinalMatrix["Motif"] = motif_label.split(".")
    FinalMatrix.set_index("Motif", inplace=True)

    # if we are in the more stringent steps of the iteration, we set flexible position rows to 0
    if iteration in ["rigid", "foldx"]:
        FinalMatrix.iloc[flex_positions] = 0

    # if we are in the iteration steps that imply using FoldX, calculate the FoldX matrix from true_positives
    # and iterate over it
    if iteration in ["foldx", "flexible"]:
        FoldXMatrix, source = foldx_process(true_positives, motif_re, motif_name, motif_label, tmp_path, pdbs_reparados_path, pdb_db_path, cwd_path)
        print(FoldXMatrix)
        print(f"Source: {source}")
        if not isinstance(FoldXMatrix, pd.DataFrame):
            return None, None, None, None, None
        if source == "pdb":
            for i in range(len(motif_label.split("."))):
                for alt in FoldXMatrix:
                    # Pathogenics, upper cutoff was determined to be 1.5
                    if FoldXMatrix.iloc[i][alt] >= 2.35:
                        if iteration == "foldx" and i not in flex_positions:
                            motif_pos = motif_label.split(".")[i]
                            for aa in motif_pos:
                                FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                        elif iteration == "flexible":
                            if i not in flex_positions:
                                motif_pos = motif_label.split(".")[i]
                                for aa in motif_pos:
                                    FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                            elif i in flex_positions:
                                for aa in aas:
                                    FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                    # Benigns, lower cutoff was determined to be 1
                    elif FoldXMatrix.iloc[i][alt] <= 1.85:
                        if iteration == "foldx" and i not in flex_positions:
                            motif_pos = motif_label.split(".")[i]
                            for aa in motif_pos:
                                FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                        elif iteration == "flexible":
                            if i not in flex_positions:
                                motif_pos = motif_label.split(".")[i]
                                for aa in motif_pos:
                                    FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                            elif i in flex_positions:
                                for aa in aas:
                                    FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                    # if we are in the final iteration step and FoldX verdict was not definitive, we apply
                    # the mean by row and see if the variant value is greater or less than it
                    elif 1.85 < FoldXMatrix.iloc[i][alt] < 2.35 and iteration == "flexible":
                        # Pathogenics
                        if FoldXMatrix.iloc[i][alt] >= FoldXMatrix.mean(axis=1).iloc[i]:
                            if i not in flex_positions:
                                motif_pos = motif_label.split(".")[i]
                                for aa in motif_pos:
                                    FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                            elif i in flex_positions:
                                for aa in aas:
                                    FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                        # Benigns
                        elif FoldXMatrix.iloc[i][alt] < FoldXMatrix.mean(axis=1).iloc[i]:
                            if i not in flex_positions:
                                motif_pos = motif_label.split(".")[i]
                                for aa in motif_pos:
                                    FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                            elif i in flex_positions:
                                for aa in aas:
                                    FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
        elif source == "af2":
            for i in range(len(motif_label.split("."))):
                for alt in FoldXMatrix:
                    # Pathogenics, upper cutoff was determined to be 1.6
                    if FoldXMatrix.iloc[i][alt] >= 1.85:
                        if iteration == "foldx" and i not in flex_positions:
                            motif_pos = motif_label.split(".")[i]
                            for aa in motif_pos:
                                FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                        elif iteration == "flexible":
                            if i not in flex_positions:
                                motif_pos = motif_label.split(".")[i]
                                for aa in motif_pos:
                                    FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                            elif i in flex_positions:
                                for aa in aas:
                                    FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                    # Benigns, lower cutoff was determined to be 1
                    elif FoldXMatrix.iloc[i][alt] <= 1.35:
                        if iteration == "foldx" and i not in flex_positions:
                            motif_pos = motif_label.split(".")[i]
                            for aa in motif_pos:
                                FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                        elif iteration == "flexible":
                            if i not in flex_positions:
                                motif_pos = motif_label.split(".")[i]
                                for aa in motif_pos:
                                    FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                            elif i in flex_positions:
                                for aa in aas:
                                    FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                    # if we are in the final iteration step and FoldX verdict was not definitive, we apply
                    # the mean by row and see if the variant value is greater or less than it
                    elif 1.35 < FoldXMatrix.iloc[i][alt] < 1.85 and iteration == "flexible":
                        # Pathogenics
                        if FoldXMatrix.iloc[i][alt] >= FoldXMatrix.mean(axis=1).iloc[i]:
                            if i not in flex_positions:
                                motif_pos = motif_label.split(".")[i]
                                for aa in motif_pos:
                                    FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                            elif i in flex_positions:
                                for aa in aas:
                                    FoldX_P.append(iteration_object_variant(aa, i, alt, "Pathogenic"))
                        # Benigns
                        elif FoldXMatrix.iloc[i][alt] < FoldXMatrix.mean(axis=1).iloc[i]:
                            if i not in flex_positions:
                                motif_pos = motif_label.split(".")[i]
                                for aa in motif_pos:
                                    FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))
                            elif i in flex_positions:
                                for aa in aas:
                                    FoldX_B.append(iteration_object_variant(aa, i, alt, "Benign"))

        # if we are in the foldx iteration step, we prioritize ClinVar and GnomAD verdicts but then we remove
        # all foldx verdicts on flexible positions, keep only foldx verdicts on rigid positions
        if iteration == "foldx":
            if source == "pdb":
                FoldXMatrix = (5/(1+(math.e**(math.log(pow(2/3, 4))*(FoldXMatrix-2.1))))) - 2.5
            elif source == "af2":
                FoldXMatrix = (5/(1+(math.e**(math.log(pow(2/3, 4))*(FoldXMatrix-1.6))))) - 2.5
            FoldXMatrix.reset_index(inplace=True, drop=True)
            FinalMatrix.reset_index(inplace=True, drop=True)
            FinalMatrix = combine_dfs([ClinVarMatrix, GnomADMatrix, FoldXMatrix], aas)
            FinalMatrix["Motif"] = motif_label.split(".")
            FinalMatrix.set_index("Motif", inplace=True)
            FinalMatrix.iloc[flex_positions] = 0

        # if we are in the flexible iteration step, we just prioritize ClinVar and GnomAD verdicts
        # and then we keep FoldX verdicts in all the other positions
        elif iteration == "flexible":
            if source == "pdb":
                FoldXMatrix = (5/(1+(math.e**(math.log(pow(2/3, 4))*(FoldXMatrix-2.1))))) - 2.5
            elif source == "af2":
                FoldXMatrix = (5/(1+(math.e**(math.log(pow(2/3, 4))*(FoldXMatrix-1.6))))) - 2.5
            FoldXMatrix.reset_index(inplace=True, drop=True)
            FinalMatrix.reset_index(inplace=True, drop=True)
            FinalMatrix = combine_dfs([ClinVarMatrix, GnomADMatrix, FoldXMatrix], aas)
            FinalMatrix["Motif"] = motif_label.split(".")
            FinalMatrix.set_index("Motif", inplace=True)
            
    else:
        source = None

    FinalMatrix = round(FinalMatrix, 1).replace(-0, 0)
    print(FinalMatrix)

    # save the FinalMatrix on tsv format
    FinalMatrix.to_csv(f"{cwd_path}Motifs/{motif_name}/FinalMatrix_{motif_name}_{iteration}_{number}.tsv", sep="\t")

    # plot the matrix
    plt.figure(figsize=(15, 5))
    cmap = sns.diverging_palette(145, 10, l=50, center="light", as_cmap=True)
    ax = sns.heatmap(FinalMatrix, vmin=-5, vmax=5, cmap=cmap, linewidths=.5, alpha=0.7, annot=True)
    sns.set(font_scale=1.2)
    plt.yticks(rotation=0)
    plt.ylabel("")
    plt.savefig(f"{cwd_path}Motifs/{motif_name}/FinalMatrix_{motif_name}_{iteration}_{number}.png")
    plt.savefig(f"{cwd_path}Motifs/{motif_name}/FinalMatrix_{motif_name}.png")
    plt.title(f"FinalMatrix, Iteration: {iteration}, Number: {number}")
    if screen == True:
        plt.show()
    else:
        plt.close()

    # priorization of ClinVar and GnomAD variants over FoldX information
    P_Variants = unique(ClinVar_P+FoldX_P)
    B_Variants = unique(ClinVar_B+GnomAD_B+FoldX_B)

#     for variant in P_Variants:
#         print(variant.__dict__)
#     for variant in B_Variants:
#         print(variant.__dict__)

    return P_Variants, B_Variants, FinalMatrix, source, true_positives

def iteration_filter(motif_name, motif_label, motif_re, iteration, screening_proteoma, true_positives, positives, negatives, remaining, Conservation_cutoff, SecStr_dict_cutoff, ExpBur_dict_cutoff, MinPos_cutoff, MaxPos_cutoff, top_GOTerms, clinvar_db_path, gnomad_db_path, uniref90_db_path, proteoma_list_path, sp_list_path, sp_fasta_path, trembl_fasta_path, tmp_path, pdbs_reparados_path, pdb_db_path, xml_db_path, sec_str_path, GO_db_path, sasa_thr_path, alphafold_human_path, scratch_path, cwd_path, plotting_dict, screen = False):
    """
    This function commands each iteration step, taking as input the remaining, positives and negatives matches
    from the previous iteration step and making the corresponding decisions to reorganize the MotSASi universe
    in new remaining, positives and negatives spaces
    """
    # variable to account if changes in the remaining group were registered, starts as sring, continues as int
    Remaining_Difference = "Starting iteration"
    # counter of repetitions in this round of the iteration step
    number = 0

    # keep on iterating while the number of remaining members change
    while Remaining_Difference != 0:

        print("*"*50)

        print(f"Iteration: {iteration}, Number: {str(number)}")

        print(f"Positives: {str(len(positives))}")
        print(f"Remaining: {str(len(remaining))}")
        print(f"Negatives: {str(len(negatives))}")
        print("Positivos, negativos, remaining")
        print(len([match for match in positives if match.ELM_Positive == True]), 
              len([match for match in negatives if match.ELM_Positive == True]),
              len([match for match in remaining if match.ELM_Positive == True]))

        # find ClinVar and GnomAD variants in true positives and positives previously found, mapped to the
        # motif positions
        P_Variants, B_Variants, FinalMatrix, source, true_positives = P_B(motif_label, motif_name, motif_re, true_positives, positives, clinvar_db_path, gnomad_db_path, cwd_path, tmp_path, pdbs_reparados_path, pdb_db_path, iteration, number, screen = screen)

        if not isinstance(FinalMatrix, pd.DataFrame):
            return None, None, None, None, None, None, None
        
        print(f"P_Variants: {str(len(P_Variants))}")
        print(f"B_Variants: {str(len(B_Variants))}")

        # list with matches whose ClinVar variants agree with ClinVar variants found in the true and new positives
        ClinVarPass = []

        # iterate through remaining and keep matches with variants that agree in the corresponding positions
        # (both benigns and pathogenics)
        for match in remaining:
            for motif_variant in [variant for variant in match.Variants if variant.source_db == "ClinVar" and variant.verdict == "Pathogenic"]:
                for P_variant in P_Variants:
                    if motif_variant.ref == P_variant.ref and \
                    motif_variant.index == P_variant.index and \
                    motif_variant.alt == P_variant.alt and \
                    motif_variant.verdict == P_variant.verdict:
                        ClinVarPass.append(match)
            for motif_variant in [variant for variant in match.Variants if variant.source_db == "ClinVar" and variant.verdict == "Benign"]:
                for B_variant in B_Variants:
                    if motif_variant.ref == B_variant.ref and \
                    motif_variant.index == B_variant.index and \
                    motif_variant.alt == B_variant.alt and \
                    motif_variant.verdict == B_variant.verdict:
                        ClinVarPass.append(match)
        print(f"ClinVarPass: {len(ClinVarPass)}")

        # list with matches whose GnomAD variants agree with GnomAD variants found in the true and new positives
        GnomADPass = []

        # iterate through remaining and keep matches with variants that agree in the corresponding positions
        # (benigns only)
        for match in remaining:
            for motif_variant in [variant for variant in match.Variants if variant.source_db == "GnomAD" and variant.verdict == "Benign"]:
                for B_variant in B_Variants:
                    if motif_variant.ref == B_variant.ref and \
                    motif_variant.index == B_variant.index and \
                    motif_variant.alt == B_variant.alt and \
                    motif_variant.verdict == B_variant.verdict:
                        GnomADPass.append(match)
        print(f"GnomADPass: {len(GnomADPass)}")
        print([match.UniProtID for match in GnomADPass])

        # sum
        TotalPass = ClinVarPass + GnomADPass
        print(f"TotalPass: {len(TotalPass)}")

        # list with matches whose ClinVar variants disagree with ClinVar variants found in the true and new positives
        ClinVarDiscard = []

        # iterate through remaining and keep matches with variants that disagree in the corresponding positions
        # (both benigns and pathogenics)
        for match in remaining:
            for motif_variant in [variant for variant in match.Variants if variant.source_db == "ClinVar" and variant.verdict == "Benign"]:
                for P_variant in P_Variants:
                    if motif_variant.ref == P_variant.ref and \
                    motif_variant.index == P_variant.index and \
                    motif_variant.alt == P_variant.alt and \
                    motif_variant.verdict != P_variant.verdict:
                        ClinVarDiscard.append(match)
            for motif_variant in [variant for variant in match.Variants if variant.source_db == "ClinVar" and variant.verdict == "Pathogenic"]:
                for B_variant in B_Variants:
                    if motif_variant.ref == B_variant.ref and \
                    motif_variant.index == B_variant.index and \
                    motif_variant.alt == B_variant.alt and \
                    motif_variant.verdict != B_variant.verdict:
                        ClinVarDiscard.append(match)
        print(f"ClinVarDiscard: {len(ClinVarDiscard)}")

        # list with matches whose GnomAD variants disagree with GnomAD variants found in the true and new positives
        GnomADDiscard = []

        # iterate through remaining and keep matches with variants that disagree in the corresponding positions
        # (benigns only)
        for match in remaining:
            for motif_variant in [variant for variant in match.Variants if variant.source_db == "GnomAD" and variant.verdict == "Benign"]:
                for P_variant in P_Variants:
                    if motif_variant.ref == P_variant.ref and \
                    motif_variant.index == P_variant.index and \
                    motif_variant.alt == P_variant.alt and \
                    motif_variant.verdict != P_variant.verdict:
                        GnomADDiscard.append(match)
        print(f"GnomADDiscard: {len(GnomADDiscard)}")

        # sum
        TotalDiscard = ClinVarDiscard + GnomADDiscard
        print(f"TotalDiscard: {len(TotalDiscard)}")

        # list with matches that present variants that both agree and disagree with the true and new positives variants
        TotalConflict = []
        for Pass in TotalPass:
            if Pass in TotalDiscard:
                if Pass not in TotalConflict:
                    Pass.conflict_label(True)
                    TotalConflict.append(Pass)
        print(f"TotalConflict: {len(TotalConflict)}")

        # for conflict motifs, we count agreements and disagreements and evaluate which ones are majority
        # and remove them from the list where they are minority
        for Conflict in TotalConflict:
#             print(Conflict.__dict__)
#             for variant in Conflict.Variants:
#                 print(variant.__dict__)
            agreements = 0
            disagreements = 0
            for Pass in TotalPass:
                if Conflict == Pass:
                    agreements += 1
            for Discard in TotalDiscard:
                if Conflict == Discard:
                    disagreements += 1
#             print("agreements, disagreements")
#             print(agreements, disagreements)
            if agreements > disagreements:
                TotalDiscard = [Discard for Discard in TotalDiscard if Discard != Conflict]
            elif agreements <= disagreements:
                TotalPass = [Pass for Pass in TotalPass if Pass != Conflict]

        # remove duplicates
        TotalPass = unique(TotalPass)
        TotalDiscard = unique(TotalDiscard)

        # label the matches of TotalPass and TotalDiscard as MotSASi positives or negatives, respectively
        for match in remaining:
            if match in TotalPass:
                match.matrix_label(True)
                match.call_MotSASi_positive()
            if match in TotalDiscard:
                match.matrix_label(False)
                match.call_MotSASi_negative()

        # candidates are those matches that passed the variant step and must be evaluated by conservation,
        # secondary structure, relative position and GO terms
        candidates = [match for match in remaining if match.MotSASi_Positive == True]
        print(f"Candidates: {str(len(candidates))}")

        # we iterate through candidates and evaluate them for each filter
        with Pool() as pool:
            output_candidates = pool.starmap(motsasi_filters, zip(candidates, repeat(motif_label), repeat(motif_name), repeat(uniref90_db_path), repeat(sp_list_path), repeat(sp_fasta_path), repeat(trembl_fasta_path), repeat(xml_db_path), repeat(sec_str_path), repeat(alphafold_human_path), repeat(scratch_path), repeat(sasa_thr_path), repeat(GO_db_path), repeat(Conservation_cutoff), repeat(SecStr_dict_cutoff), repeat(ExpBur_dict_cutoff), repeat(MinPos_cutoff), repeat(MaxPos_cutoff), repeat(top_GOTerms), repeat(cwd_path), repeat(tmp_path)))

        for i in range(len(remaining)):
            for candidate in output_candidates:
                if remaining[i].UniProtID == candidate.UniProtID and remaining[i].start == candidate.start:
                    #print(candidate.__dict__)
                    remaining[i] = candidate
        for i in range(len(true_positives)):
            for candidate in output_candidates:
                if true_positives[i].UniProtID == candidate.UniProtID and true_positives[i].start == candidate.start:
                    #print(candidate.__dict__)
                    true_positives[i] = candidate

        # matchs that passed all filters
        FinalCand_Number = len([match for match in remaining if match.MotSASi_Positive == True and \
                    match.passed_conservation == True and \
                    match.passed_secondary_structure == True and \
                    match.passed_relative_position == True and \
                    match.passed_GOTerms == True and \
                    match.passed_exposure == True])

        # we include in the positives group all new positives
        positives += [match for match in remaining if match.MotSASi_Positive == True and \
                    match.passed_conservation == True and \
                    match.passed_secondary_structure == True and \
                    match.passed_relative_position == True and \
                    match.passed_GOTerms == True and \
                    match.passed_exposure == True]

        # we revert the MotSASi label in those matches that passed the variant filter but then did not make it
        for match in remaining:
            if match.MotSASi_Positive == True and \
                    (match.passed_conservation == False or \
                    match.passed_secondary_structure == False or \
                    match.passed_relative_position == False or \
                    match.passed_GOTerms == False or \
                    match.passed_exposure == False):
                match.revert_MotSASi_verdict()

        # we include in the negatives group all new negatives
        negatives += [match for match in remaining if match.MotSASi_Negative == True]

        # we calculate the initial amount of matches in remaining
        Initial_Rem_Length = len(remaining)
        print(f"Longitud inicial de remaining: {Initial_Rem_Length}")
        # we update the remaining group
        remaining = [match for match in remaining if match.MotSASi_Positive == None and match.MotSASi_Negative == None]
        # we calculate the final amount of matches in remaining and make the difference
        Final_Rem_Length = len(remaining)
        print(f"Longitud final de remaining: {Final_Rem_Length}")
        Remaining_Difference = Initial_Rem_Length - Final_Rem_Length

        print(f"Positives: {str(len(positives))}")
        print(f"Negatives: {str(len(negatives))}")

        # extension of the MotSASi universe (should be constant)
        MotSASi_Universe = len(positives) + len(negatives) + len(remaining)

        # matchs of the MotSASi universe that were accepted or rejected
        MotSASi_Universe_resolved = len(positives) + len(negatives)

        # candidates that were evaluated in this round
        Cand_Number = len(candidates)

        # this is the amount of true positives that were relevated by MotSASi and remain in any group
        Total_ELM_Positives = len(true_positives)
        Total_ELM_Positives_omim = len([match for match in true_positives if match.OMIM == True])

        # true positives that were relevated by MotSASi and could be evaluated until this moment
        ELM_Pos_MotSASi = len([match for match in positives if match.ELM_Positive == True]) + \
                            len([match for match in negatives if match.ELM_Positive == True])

        # true positives that were relevated by MotSASi and could be evaluated until this moment and passed
        # all filters
        ELM_Pos_MotSASi_Pos = len([match for match in positives if match.ELM_Positive == True])

        # true positives that were relevated by MotSASi and could be evaluated until this moment and did not
        # pass all filters
        ELM_Pos_MotSASi_Neg = len([match for match in negatives if match.ELM_Positive == True])
        
        # true positives that were relevated by MotSASi and could be evaluated until this moment and passed
        # all filters but with conflicts
        ELM_Pos_MotSASi_Pos_conflict = len([match for match in positives if match.ELM_Positive == True and match.conflict == True])

        # true positives that were relevated by MotSASi and could be evaluated until this moment
        ELM_Pos_with_bp_variants = len([match for match in positives if match.ELM_Positive == True]) + \
                            len([match for match in negatives if match.ELM_Positive == True]) + \
                            len([match for match in remaining if match.ELM_Positive == True])

        # matches that were discarded by the MotSASi matrix
        Matrix_Discarded = len([match for match in negatives if match.passed_matrix == False])

        # matches that passed the MotSASi matrix
        Matrix_Passed = len([match for match in negatives if match.passed_matrix == True]) + \
                        len([match for match in positives if match.passed_matrix == True])

        # matches that suffered from conflicts in the matrix
        Matrix_conflicts = len([match for match in negatives if match.conflict == True]) + \
                        len([match for match in positives if match.conflict == True])

        # matches that contain ClinVar variants
        clinvar_variants_matches = [match for match in screening_proteoma if len([variant for variant in match.Variants if variant.source_db == "ClinVar"]) > 0]

        # matches that contain ClinVar BP variants
        clinvar_bp_variants_matches = [match for match in screening_proteoma if len([variant for variant in match.Variants if variant.source_db == "ClinVar" and variant.verdict in ["Benign", "Pathogenic"]]) > 0]

        # matches that contain ClinVar variants
        gnomad_variants_matches = [match for match in screening_proteoma if len([variant for variant in match.Variants if variant.source_db == "GnomAD"]) > 0]

        # matches that contain ClinVar BP variants
        gnomad_b_variants_matches = [match for match in screening_proteoma if len([variant for variant in match.Variants if variant.source_db == "GnomAD" and variant.verdict in ["Benign"]]) > 0]

#         if len(Missings_Cons) > 0:
#             print("Some Conservation scores have no values... check Missings_Cons.txt file")
#         if len(Missings_SecStr) > 0:
#             print("Some Secondary Structures have no values... check Missings_SecStr.txt file")
#         if len(Missings_GOTerms) > 0:
#             print("Some GO Terms have no values... check Missings_GOTerms.txt file")
#         if len(Missings_ExpBur) > 0:
#             print("Some Exposures have no values... check Missings_ExpBur.txt file")

        print(f"Results of Iteration: {iteration}, Number: {number}")
        print(f"MotSASi Universe: {MotSASi_Universe}")
        print(f"MotSASi Universe with resolution: {MotSASi_Universe_resolved}")
        print(f"Candidates in this round (passed the SAS matrix): {Cand_Number}")
        print(f"Final Candidates: {FinalCand_Number}")
        print(f"ELM Positives: {Total_ELM_Positives}")
        print(f"ELM Positives analyzables by MotSASi: {ELM_Pos_MotSASi}")
        print(f"ELM Positives detected: {ELM_Pos_MotSASi_Pos}")
        print("Positivos, negativos, remaining")
        print(len([match for match in positives if match.ELM_Positive == True]), 
              len([match for match in negatives if match.ELM_Positive == True]),
              len([match for match in remaining if match.ELM_Positive == True]))
        
        
        if Remaining_Difference > 0 or number == 0:
            # plotting purpose
            iter_plotting_dict = {}
            iter_plotting_dict["MotSASi_ELM_Positives"] = ELM_Pos_with_bp_variants
            positives_less_ELM = len([match for match in positives if match.ELM_Positive == None])
            iter_plotting_dict["MotSASi_Positives"] = positives_less_ELM
            negatives_less_ELM = len([match for match in negatives if match.ELM_Positive == None])
            iter_plotting_dict["MotSASi_Negatives"] = negatives_less_ELM
            remaining_less_ELM = len([match for match in remaining if match.ELM_Positive == None])
            iter_plotting_dict["MotSASi_Remaining"] = remaining_less_ELM
            no_actionable_variants_motifs = len([motif for motif in screening_proteoma if len(motif.Variants) > 0 and motif.number_MotSASi_variants == 0])
            no_actionable_variants_motifs_noELM = len([motif for motif in screening_proteoma if len(motif.Variants) > 0 and motif.number_MotSASi_variants == 0 and motif.ELM_Positive == None])
            no_actionable_variants_motifs_noELM_omim = len([motif for motif in screening_proteoma if len(motif.Variants) > 0 and motif.number_MotSASi_variants == 0 and motif.OMIM == True and motif.ELM_Positive == None])
            iter_plotting_dict["No_Actionable_Variants_Motifs"] = no_actionable_variants_motifs
            no_variants_motifs = len([motif for motif in screening_proteoma if len(motif.Variants) == 0])
            no_variants_motifs_noELM = len([motif for motif in screening_proteoma if len(motif.Variants) == 0 and motif.ELM_Positive == None])
            no_variants_motifs_noELM_omim = len([motif for motif in screening_proteoma if len(motif.Variants) == 0 and motif.OMIM == True and motif.ELM_Positive == None])
            iter_plotting_dict["No_Variants_Motifs"] = no_variants_motifs
            plotting_dict[f"{iteration}_{number}"] = iter_plotting_dict

        # number of the next round inside this iteration step
        number += 1

    motif_statistics = {"Motif": motif_name, "True_Positives": ELM_Pos_MotSASi_Pos, "False_Negatives": ELM_Pos_MotSASi_Neg,
                        "Conflicts_Positives": ELM_Pos_MotSASi_Pos_conflict, "Total_Positives_with_bp_variants": ELM_Pos_with_bp_variants, 
                        "Total_ELM_Positives": Total_ELM_Positives, "Total_ELM_Positives_OMIM": Total_ELM_Positives_omim,
                        "Discarded_with_the_matrix": Matrix_Discarded, 
                        "Pass_the_matrix": Matrix_Passed, "Conflicts": Matrix_conflicts,
                        "Pos_neg_motifs_with_variants": MotSASi_Universe_resolved, "Motifs_with_variants": MotSASi_Universe, 
                        "Total_hits_in_the_proteome": len(screening_proteoma),
                        "Total_hits_in_the_proteome_noELM": len([match for match in screening_proteoma if match.ELM_Positive == None]),
                        "Total_hits_in_the_proteome_noELM_OMIM": len([match for match in screening_proteoma if match.ELM_Positive == None and match.OMIM == True]),
                        "Sensitivity": ELM_Pos_MotSASi_Pos/Total_ELM_Positives, 
                        "To_Analyze_with_MotSASi": (MotSASi_Universe-Total_ELM_Positives)/MotSASi_Universe, 
                        "Initially_Discarded_in_MotSASi": Matrix_Discarded/MotSASi_Universe, 
                        "To_Analyze_with_Sequence_based": Matrix_Passed/MotSASi_Universe,
                        "MotSASi_Positives": positives_less_ELM, "MotSASi_Negatives": negatives_less_ELM,
                        "MotSASi_Remaining": remaining_less_ELM, 
                        "No_Actionable_Variants_Motifs": no_actionable_variants_motifs,
                        "No_Variants_Motifs": no_variants_motifs,
                        "No_Actionable_Variants_Motifs_noELM": no_actionable_variants_motifs_noELM,
                        "No_Variants_Motifs_noELM": no_variants_motifs_noELM,
                        "No_Actionable_Variants_Motifs_noELM_OMIM": no_actionable_variants_motifs_noELM_omim,
                        "No_Variants_Motifs_noELM_OMIM": no_variants_motifs_noELM_omim,
                        "ClinVar_Variants_Matches": len(clinvar_variants_matches),
                        "ClinVar_bp_Variants_Matches": len(clinvar_bp_variants_matches),
                        "GnomAD_Variants_Matches": len(gnomad_variants_matches),
                        "GnomAD_b_Variants_Matches": len(gnomad_b_variants_matches),
                        "Source": source}

    return true_positives, positives, negatives, remaining, plotting_dict, motif_statistics, FinalMatrix

def motsasi_filters(candidate, motif_label, motif_name, uniref90_db_path, sp_list_path, sp_fasta_path, trembl_fasta_path, xml_db_path, sec_str_path, alphafold_human_path, scratch_path, sasa_thr_path, GO_db_path, Conservation_cutoff, SecStr_dict_cutoff, ExpBur_dict_cutoff, MinPos_cutoff, MaxPos_cutoff, top_GOTerms, cwd_path, tmp_path):

    # Conservation filter
    try:
        Cand_Cons = ConservationScore([candidate], motif_label, motif_name, uniref90_db_path, sp_list_path, sp_fasta_path, trembl_fasta_path, cwd_path, tmp_path, plot = False, save_tsv = False)
        if Cand_Cons[candidate.UniProtID].mean() >= Conservation_cutoff:
            candidate.conservation_label(True)
        elif Cand_Cons[candidate.UniProtID].mean() < Conservation_cutoff:
            candidate.conservation_label(False)
    except:
        print(f"Falla la conservación: {candidate.UniProtID}")
        candidate.conservation_label(False)
        with open(f"{cwd_path}Motifs/{motif_name}/Missings_Cons.txt", "a") as output:
            output.write(f"{candidate.UniProtID}\n")
            output.close()

    # Secondary Structure filter
    try:
        Cand_SecStr = scratch([candidate], motif_name, motif_label, xml_db_path, sec_str_path, alphafold_human_path, scratch_path, cwd_path, tmp_path, plot = False, save_tsv = False, pred_number = 3)
        agreements = 0
        disagreements = 0
        for n, (i, row) in enumerate(Cand_SecStr.iterrows()):
            for secstr in Cand_SecStr:
                if Cand_SecStr.iloc[n][secstr] == 100:
                    if secstr in SecStr_dict_cutoff[n]:
                        agreements += 1
                    elif secstr not in SecStr_dict_cutoff[n]:
                        disagreements += 1
        if agreements >= disagreements:
            candidate.secondary_structure_label(True)
        elif agreements < disagreements:
            candidate.secondary_structure_label(False)
    except:
        print(f"Falla la estructura secundaria: {candidate.UniProtID}")
        candidate.secondary_structure_label(False)
        with open(f"{cwd_path}Motifs/{motif_name}/Missings_SecStr.txt", "a") as output:
            output.write(f"{candidate.UniProtID}\n")
            output.close()

    # Exposure filter
    try:
        Cand_ExpBur = sasa([candidate], motif_name, motif_label, xml_db_path, alphafold_human_path, sasa_thr_path, cwd_path, tmp_path, plot = False, save_tsv = False)
        agreements = 0
        disagreements = 0
        for n, (i, row) in enumerate(Cand_ExpBur.iterrows()):
            for expbur in Cand_ExpBur:
                if Cand_ExpBur.iloc[n][expbur] == 100:
                    if expbur in ExpBur_dict_cutoff[n]:
                        agreements += 1
                    elif expbur not in ExpBur_dict_cutoff[n]:
                        disagreements += 1
        if agreements >= disagreements:
            candidate.exposure_label(True)
        elif agreements < disagreements:
            candidate.exposure_label(False)
    except:
        print(f"Falla la exposición: {candidate.UniProtID}")
        candidate.exposure_label(False)
        with open(f"{cwd_path}Motifs/{motif_name}/Missings_ExpBur.txt", "a") as output:
            output.write(f"{candidate.UniProtID}\n")
            output.close()

    # Relative Position filter
#     Cand_RelPos = ((candidate.start + candidate.end)//2)/candidate.length
#     if MinPos_cutoff <= Cand_RelPos <= MaxPos_cutoff:
    candidate.relative_position_label(True)
#     else:
#         candidate.relative_position_label(False)

    # GO Terms filter
    try:
        Cand_GOTerms = Motif_GO([candidate], motif_name, xml_db_path, GO_db_path, cwd_path, plot=False, save_tsv=False).columns.to_list()
        if len([GOTerm for GOTerm in Cand_GOTerms if GOTerm in top_GOTerms]) >= 1:
            candidate.GOTerms_label(True)
        elif len([GOTerm for GOTerm in Cand_GOTerms if GOTerm in top_GOTerms]) == 0:
            candidate.GOTerms_label(False)
    except:
        print(f"Fallan los GO: {candidate.UniProtID}")
        candidate.GOTerms_label(False)
        with open(f"{cwd_path}Motifs/{motif_name}/Missings_GOTerms.txt", "a") as output:
            output.write(f"{candidate.UniProtID}\n")
            output.close()

    return candidate

def ConsCutOff(true_positives, motif_label, motif_name, uniref90_db_path, proteoma_list_path, sp_list_path, sp_fasta_path, trembl_fasta_path, cwd_path, tmp_path, screen):
    """
    This function takes the true positives as input and outputs the minimum conservation score observed
    along the motif
    """
    # dataframe with motif residues in rows and true positives in columns
    Conservation_df = ConservationScore(true_positives, motif_label, motif_name, uniref90_db_path, sp_list_path, sp_fasta_path, trembl_fasta_path, cwd_path, tmp_path, screen = screen)
    # we calculate the mean by row (by motif residue) and then take the minimum of them as cut off
    Conservation_cutoff = Conservation_df.mean(axis=1, numeric_only = True).min()

    return Conservation_cutoff

def SecStrCutOff(true_positives, motif_name, motif_label, xml_db_path, sec_str_path, alphafold_human_path, scratch_path, cwd_path, tmp_path, screen):
    """
    This function takes the true positives as input and outputs a dictionary where keys are motif residues
    and values are secondary structure elements allowed in each position
    """
    # dataframe with motif residues in rows and secondary structure elements seen in true positives
    # in columns (as percentage of ocurrences)
    SecStr_df = scratch(true_positives, motif_name, motif_label, xml_db_path, sec_str_path, alphafold_human_path, scratch_path, cwd_path, tmp_path, pred_number = 3, screen = screen)

    # dictionary with elements(future output)
    SecStr_dict_cutoff = {}
    
    # we iterate through rows and keep elements that appear al least once in a true positive
    for n, (i, row) in enumerate(SecStr_df.iterrows()):
        potential = []
        for secstr in SecStr_df:
            if SecStr_df.iloc[n][secstr] > 5:
                potential.append(secstr)
        SecStr_dict_cutoff[n] = potential
    
    return SecStr_dict_cutoff

def ExpBurCutOff(true_positives, motif_name, motif_label, xml_db_path, alphafold_human_path, sasa_thr_path, cwd_path, tmp_path, screen):
    """
    This function takes the true positives as input and outputs a dictionary where keys are motif residues
    and values are exposures allowed in each position
    """
    # dataframe with motif residues in rows and exposures seen in true positives in alphafold
    # in columns (as percentage of ocurrences)
    ExpBur_df = sasa(true_positives, motif_name, motif_label, xml_db_path, alphafold_human_path, sasa_thr_path, cwd_path, tmp_path, screen = screen)

    # dictionary with elements(future output)
    ExpBur_dict_cutoff = {}

    # we iterate through rows and keep elements that appear al least once in a true positive
    for n, (i, row) in enumerate(ExpBur_df.iterrows()):
        potential = []
        for expbur in ExpBur_df:
            if ExpBur_df.iloc[n][expbur] > 5:
                potential.append(expbur)
        ExpBur_dict_cutoff[n] = potential
    
    return ExpBur_dict_cutoff

def GOTermsCutOff(true_positives, motif_name, xml_db_path, GO_db_path, cwd_path, screen):
    """
    This function takes the true positives as input and outputs a list with a certain amount of representative GO terms
    """
    # dataframe with GO terms in columns and percentage of ocurrences in the only row
    GO_df = Motif_GO(true_positives, motif_name, xml_db_path, GO_db_path, cwd_path, screen = screen)
    GO_df = GO_df.T
    # try to keep very representative terms, if not enough take a minimum of 5 terms
    top_GOTerms = GO_df[GO_df["%"] >= 50].index.tolist()
    if len(top_GOTerms) < 5:
        GO_df = GO_df[GO_df["%"] >= 10]
        top_GOTerms = GO_df.iloc[:5].index.tolist()

    return top_GOTerms

def unique(objects_list):
    """
    This function has the only purpose of removing duplicates from objects lists
    """
    output_list = []
    for ind_object in objects_list:
        if ind_object not in output_list:
            output_list.append(ind_object)
    return output_list

def iteration_process_plot(plotting_dict, motif_name, cwd_path, screen=False):
    
    # Convertir el diccionario a un DataFrame
    df = pd.DataFrame(plotting_dict)

    # Configurar el ancho de las barras
    ancho_barra = 0.5

    # Inicializar los valores acumulados (para la pila)
    valores_acumulados = np.zeros(len(df.columns))

    # Crear el gráfico
    fig = plt.figure(figsize = (12, 7))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor("lightseagreen")
    plt.gcf().set_facecolor("lightgrey")

    # Colores para cada subcategoría
    colores = ["violet", "limegreen", "red", "grey", "lightblue", "royalblue"]

    # Dibujar cada fila (subcategoría)
    for i, fila in enumerate(df.index):
        ax.bar(df.columns, df.loc[fila], ancho_barra, bottom=valores_acumulados, label=fila, color=colores[i])

        # Actualizar los valores acumulados para la próxima pila
        valores_acumulados += df.loc[fila].values

    # Añadir etiquetas y título
    ax.set_ylabel("")
    plt.xticks(rotation="horizontal")
    ax.set_title("Iteration process progression")
    ax.legend(title="Categories", bbox_to_anchor=(1.0, 1), loc="upper left")

    plt.savefig(f"{cwd_path}Motifs/{motif_name}/Iteration_process_plot_{motif_name}.png", bbox_inches = "tight")

    # Show
    if screen == True:
        plt.show()
    else:
        plt.close()

    return

def variants_prediction(motif_name, true_positives, positives, FinalMatrix, cwd_path):

    if f"motsasi_variant_predictions.csv" not in os.listdir(cwd_path):
        df = pd.DataFrame(columns =["Motif_name", "UniProtID", "Start", "Ref", "Pos", "Alt", "MotSASi_score", "ELM_Positive", "OMIM_Phenotype"])
        df.to_csv(f"{cwd_path}motsasi_variant_predictions.csv", index = False)
        processed = []

    else:
        df = pd.read_csv(f"{cwd_path}motsasi_variant_predictions.csv")
        processed = []
        for n, row in df.iterrows():
            processed.append(f"{row.Motif_name}-{row.UniProtID}-{row.Start}-{row.Ref}-{row.Pos}-{row.Alt}-{row.MotSASi_score}-{row.ELM_Positive}-{row.OMIM_Phenotype}")
        processed = set(processed)

    FinalMatrix = FinalMatrix.reset_index(drop=True)

    true_and_new_positives = unique(true_positives + positives)

    new_variants = []

    for match in true_and_new_positives:
        for n, row in FinalMatrix.iterrows():
            for aa in FinalMatrix:
                if match.seq[n] != aa:
                    MotSASi_score = FinalMatrix.loc[n, aa]
                    if f"{motif_name}-{match.UniProtID}-{match.start}-{match.seq[n]}-{match.start+n}-{aa}-{MotSASi_score}-{match.ELM_Positive}-{match.OMIM}" not in processed:
                        if match.ELM_Positive == True:
                            aux_dict = {"Motif_name": motif_name, "UniProtID": match.UniProtID, "Start": match.start, "Ref": match.seq[n], "Pos": match.start+n,
                                       "Alt": aa, "MotSASi_score": MotSASi_score, "ELM_Positive": True, "OMIM_Phenotype": match.OMIM}
                        else:
                            aux_dict = {"Motif_name": motif_name, "UniProtID": match.UniProtID, "Start": match.start, "Ref": match.seq[n], "Pos": match.start+n,
                                       "Alt": aa, "MotSASi_score": MotSASi_score, "ELM_Positive": False, "OMIM_Phenotype": match.OMIM}
                        new_variants.append(aux_dict)

    for new_variant in new_variants:
        df = df._append(new_variant, ignore_index=True)

    df.to_csv(f"{cwd_path}motsasi_variant_predictions.csv", index=False)

    return

def negative_motifs(motif_name, negatives, cwd_path):

    if f"motsasi_negative_motifs.csv" not in os.listdir(cwd_path):
        df = pd.DataFrame(columns =["Motif_name", "UniProtID", "Start", "OMIM_Phenotype", "ELM_Positive"])
        df.to_csv(f"{cwd_path}motsasi_negative_motifs.csv", index = False)
        processed = []

    else:
        df = pd.read_csv(f"{cwd_path}motsasi_negative_motifs.csv")
        processed = []
        for n, row in df.iterrows():
            processed.append(f"{row.Motif_name}-{row.UniProtID}-{row.Start}-{row.OMIM_Phenotype}-{row.ELM_Positive}")
        processed = set(processed)

    neg_motifs = []

    for match in negatives:
        if f"{motif_name}-{match.UniProtID}-{match.start}-{match.OMIM}-{match.ELM_Positive}" not in processed:
            aux_dict = {"Motif_name": motif_name, "UniProtID": match.UniProtID, "Start": match.start, "OMIM_Phenotype": match.OMIM, "ELM_Positive": match.ELM_Positive}
            neg_motifs.append(aux_dict)

    for neg_motif in neg_motifs:
        df = df._append(neg_motif, ignore_index=True)

    df.to_csv(f"{cwd_path}motsasi_negative_motifs.csv", index=False)

    return

def iteration_process(motif_name, motif_re, motif_label, clinvar_db_path, gnomad_db_path, sp_h_genes_af_cut_off_path, uniref90_db_path, proteoma_list_path, sp_list_path, sp_fasta_path, trembl_fasta_path, tmp_path, fasta_db_path, xml_db_path, sec_str_path, GO_db_path, elm_db_path, pdb_db_path, pdbs_reparados_path, alphafold_human_path, sasa_thr_path, cwd_path, scratch_path, uniprot_omim_phen_path, screen = False):
    """
    This is the function that orders the iteration steps, running rigid, foldx, and flexible steps in order
    """
    if f"motsasi_statistics.csv" not in os.listdir(f"{cwd_path}"):
        df = pd.DataFrame(columns =["Motif", "True_Positives", "False_Negatives", "Conflicts_Positives", 
                                   "Total_Positives_with_bp_variants", "Total_ELM_Positives", "Total_ELM_Positives_OMIM",
                                   "Discarded_with_the_matrix", "Pass_the_matrix", "Conflicts", "Pos_neg_motifs_with_variants",
                                   "Motifs_with_variants", "Total_hits_in_the_proteome", "Total_hits_in_the_proteome_noELM", "Total_hits_in_the_proteome_noELM_OMIM",
                                   "Sensitivity", "To_Analyze_with_MotSASi", "Initially_Discarded_in_MotSASi",
                                   "To_Analyze_with_Sequence_based", "MotSASi_Positives", "MotSASi_Negatives", 
                                    "MotSASi_Remaining", "No_Actionable_Variants_Motifs",
                                    "No_Variants_Motifs", 
                                    "No_Actionable_Variants_Motifs_noELM", "No_Variants_Motifs_noELM",
                                    "No_Actionable_Variants_Motifs_noELM_OMIM", "No_Variants_Motifs_noELM_OMIM",
                                    "ClinVar_Variants_Matches",
                                    "ClinVar_bp_Variants_Matches",
                                    "GnomAD_Variants_Matches",
                                    "GnomAD_b_Variants_Matches", "Source"])
        df.to_csv(f"{cwd_path}motsasi_statistics.csv", index = False)
    # check if the motif has already been processed
    elif f"motsasi_statistics.csv" in os.listdir(cwd_path):
        df = pd.read_csv(f"{cwd_path}motsasi_statistics.csv", index_col=None)
        already_done = set(df.Motif.to_list())
        print(already_done)
        if motif_name in already_done:
            print(f"Motif already analyzed: {motif_name}")
            return None
        elif motif_name not in already_done:
            df = df[df.Motif != motif_name]

    # we first look for all possible matches of the motif, label the true and false ELM positives and find
    # ClinVar and GnomAD variants for them, remember that our universe will be delimited by potential motifs
    # carrying variants
    screening_proteoma = HumanProteomeSearch(motif_name, motif_re, fasta_db_path, xml_db_path, proteoma_list_path, uniprot_omim_phen_path)
    positive_ctrl(screening_proteoma, motif_name, pdb_db_path, elm_db_path)
    ClinVar_variant_objects, screening_proteoma = ClinVarMotif(screening_proteoma, motif_name, clinvar_db_path)
    GnomAD_variant_objects, screening_proteoma = GnomADMotif(screening_proteoma, motif_name, gnomad_db_path, sp_h_genes_af_cut_off_path)
    true_positives = [match for match in screening_proteoma if match.ELM_Positive == True]
    print("Motivos con variantes")
    print(len([match for match in screening_proteoma if len([variant for variant in match.Variants if variant.verdict in ["Benign", "Pathogenic"]]) > 0]))
    print(f"TRUE POSITIVES: {len(true_positives)}")
#     for motif in true_positives:
#         print(motif.__dict__)
#         for variant in motif.Variants:
#             if variant.verdict in ["Pathogenic", "Benign"]:
#                 print(variant.__dict__)

    if len(true_positives) == 0:
        print("True positives not currently available in human")
        motif_statistics = {"Motif": motif_name, "True_Positives": "No_TP", "False_Negatives": "No_TP", 
                        "Conflicts_Positives": "No_TP", "Total_Positives_with_bp_variants": "No_TP", 
                        "Total_ELM_Positives": "No_TP", "Total_ELM_Positives_OMIM": "No_TP", "Discarded_with_the_matrix": "No_TP", 
                        "Pass_the_matrix": "No_TP", "Conflicts": "No_TP",
                        "Pos_neg_motifs_with_variants": "No_TP", "Motifs_with_variants": "No_TP", 
                        "Total_hits_in_the_proteome": "No_TP",
                        "Total_hits_in_the_proteome_noELM": "No_TP",
                        "Total_hits_in_the_proteome_noELM_OMIM": "No_TP",
                        "Sensitivity": "No_TP", 
                        "To_Analyze_with_MotSASi": "No_TP", 
                        "Initially_Discarded_in_MotSASi": "No_TP", 
                        "To_Analyze_with_Sequence_based": "No_TP",
                        "MotSASi_Positives": "No_TP", "MotSASi_Negatives": "No_TP",
                        "MotSASi_Remaining": "No_TP", 
                        "No_Actionable_Variants_Motifs": "No_TP",
                        "No_Variants_Motifs": "No_TP",
                        "No_Actionable_Variants_Motifs_noELM": "No_TP", "No_Variants_Motifs_noELM": "No_TP",
                        "No_Actionable_Variants_Motifs_noELM_OMIM": "No_TP", "No_Variants_Motifs_noELM_OMIM": "No_TP",
                        "ClinVar_Variants_Matches": "No_TP",
                        "ClinVar_bp_Variants_Matches": "No_TP",
                        "GnomAD_Variants_Matches": "No_TP",
                        "GnomAD_b_Variants_Matches": "No_TP",
                        "Source": "No_TP"}
        df.loc[len(df)] = motif_statistics
        df.to_csv(f"{cwd_path}motsasi_statistics.csv", index=False)
        return None

    positives = []
    negatives = []
    #remaining = [match for match in screening_proteoma if len(match.Variants) > 0]
    remaining = [match for match in screening_proteoma if match.number_MotSASi_variants > 0]
    print(f"{str(len(positives))} motifs in P0")
    print(f"{str(len(negatives))} motifs in N0")
    print(f"{str(len(remaining))} motifs in R0")
    
    if len(remaining) == 0:
        print("There are no motif matches carrying ClinVar/gnomAD variants")
        motif_statistics = {"Motif": motif_name, "True_Positives": "No_Variants", "False_Negatives": "No_Variants", 
                        "Conflicts_Positives": "No_Variants", "Total_Positives_with_bp_variants": "No_Variants", 
                        "Total_ELM_Positives": "No_Variants", "Total_ELM_Positives_OMIM": "No_Variants", "Discarded_with_the_matrix": "No_Variants", 
                        "Pass_the_matrix": "No_Variants", "Conflicts": "No_Variants",
                        "Pos_neg_motifs_with_variants": "No_Variants", "Motifs_with_variants": "No_Variants", 
                        "Total_hits_in_the_proteome": "No_Variants",
                        "Total_hits_in_the_proteome_noELM": "No_Variants",
                        "Total_hits_in_the_proteome_noELM_OMIM": "No_Variants",
                        "Sensitivity": "No_Variants", 
                        "To_Analyze_with_MotSASi": "No_Variants", 
                        "Initially_Discarded_in_MotSASi": "No_Variants", 
                        "To_Analyze_with_Sequence_based": "No_Variants",
                        "MotSASi_Positives": "No_Variants", "MotSASi_Negatives": "No_Variants",
                        "MotSASi_Remaining": "No_Variants", 
                        "No_Actionable_Variants_Motifs": "No_Variants",
                        "No_Variants_Motifs": "No_Variants",
                        "No_Actionable_Variants_Motifs_noELM": "No_Variants", "No_Variants_Motifs_noELM": "No_Variants",
                        "No_Actionable_Variants_Motifs_noELM_OMIM": "No_Variants", "No_Variants_Motifs_noELM_OMIM": "No_Variants",
                        "ClinVar_Variants_Matches": "No_Variants",
                        "ClinVar_bp_Variants_Matches": "No_Variants",
                        "GnomAD_Variants_Matches": "No_Variants",
                        "GnomAD_b_Variants_Matches": "No_Variants",
                        "Source": "No_Variants"}
        df.loc[len(df)] = motif_statistics
        df.to_csv(f"{cwd_path}motsasi_statistics.csv", index=False)
        return None

    # plotting dictionary
    plotting_dict = {}

    # setting our conservation, secondary structure, relative position and GO terms cut offs
    Conservation_cutoff = ConsCutOff(true_positives, motif_label, motif_name, uniref90_db_path, proteoma_list_path, sp_list_path, sp_fasta_path, trembl_fasta_path, cwd_path, tmp_path, screen)
    SecStr_dict_cutoff = SecStrCutOff(true_positives, motif_name, motif_label, xml_db_path, sec_str_path, alphafold_human_path, scratch_path, cwd_path, tmp_path, screen)
    top_GOTerms = GOTermsCutOff(true_positives, motif_name, xml_db_path, GO_db_path, cwd_path, screen)
    ExpBur_dict_cutoff = ExpBurCutOff(true_positives, motif_name, motif_label, xml_db_path, alphafold_human_path, sasa_thr_path, cwd_path, tmp_path, screen)
    print(ExpBur_dict_cutoff)

    # rigid iteration
    true_positives, positives, negatives, remaining, plotting_dict, motif_statistics, FinalMatrix = iteration_filter(motif_name, motif_label, motif_re, "rigid", screening_proteoma, true_positives, positives, negatives, remaining, Conservation_cutoff, SecStr_dict_cutoff, ExpBur_dict_cutoff, MinPos_cutoff, MaxPos_cutoff, top_GOTerms, clinvar_db_path, gnomad_db_path, uniref90_db_path, proteoma_list_path, sp_list_path, sp_fasta_path, trembl_fasta_path, tmp_path, pdbs_reparados_path, pdb_db_path, xml_db_path, sec_str_path, GO_db_path, sasa_thr_path, alphafold_human_path, scratch_path, cwd_path, plotting_dict, screen = screen)
    print(f"{str(len(positives))} motifs in P1f")
    print(f"{str(len(negatives))} motifs in N1f")
    print(f"{str(len(remaining))} motifs in R1f")
    print(plotting_dict)

    # foldx iteration
    true_positives, positives, negatives, remaining, plotting_dict, motif_statistics, FinalMatrix = iteration_filter(motif_name, motif_label, motif_re, "foldx", screening_proteoma, true_positives, positives, negatives, remaining, Conservation_cutoff, SecStr_dict_cutoff, ExpBur_dict_cutoff, MinPos_cutoff, MaxPos_cutoff, top_GOTerms, clinvar_db_path, gnomad_db_path, uniref90_db_path, proteoma_list_path, sp_list_path, sp_fasta_path, trembl_fasta_path, tmp_path, pdbs_reparados_path, pdb_db_path, xml_db_path, sec_str_path, GO_db_path, sasa_thr_path, alphafold_human_path, scratch_path, cwd_path, plotting_dict, screen = screen)
    if not isinstance(FinalMatrix, pd.DataFrame):
        print("Stability matrices not available, review AF2 output")
        motif_statistics = {"Motif": motif_name, "True_Positives": "No_AF2", "False_Negatives": "No_AF2", 
                        "Conflicts_Positives": "No_AF2", "Total_Positives_with_bp_variants": "No_AF2", 
                        "Total_ELM_Positives": "No_AF2", "Total_ELM_Positives_OMIM": "No_AF2", "Discarded_with_the_matrix": "No_AF2", 
                        "Pass_the_matrix": "No_AF2", "Conflicts": "No_AF2",
                        "Pos_neg_motifs_with_variants": "No_AF2", "Motifs_with_variants": "No_AF2", 
                        "Total_hits_in_the_proteome": "No_AF2",
                        "Total_hits_in_the_proteome_noELM": "No_AF2",
                        "Total_hits_in_the_proteome_noELM_OMIM": "No_AF2",
                        "Sensitivity": "No_AF2", 
                        "To_Analyze_with_MotSASi": "No_AF2", 
                        "Initially_Discarded_in_MotSASi": "No_AF2", 
                        "To_Analyze_with_Sequence_based": "No_AF2",
                        "MotSASi_Positives": "No_AF2", "MotSASi_Negatives": "No_AF2",
                        "MotSASi_Remaining": "No_AF2", 
                        "No_Actionable_Variants_Motifs": "No_AF2",
                        "No_Variants_Motifs": "No_AF2",
                        "No_Actionable_Variants_Motifs_noELM": "No_AF2", "No_Variants_Motifs_noELM": "No_AF2",
                        "No_Actionable_Variants_Motifs_noELM_OMIM": "No_AF2", "No_Variants_Motifs_noELM_OMIM": "No_AF2",
                        "ClinVar_Variants_Matches": "No_AF2",
                        "ClinVar_bp_Variants_Matches": "No_AF2",
                        "GnomAD_Variants_Matches": "No_AF2",
                        "GnomAD_b_Variants_Matches": "No_AF2",
                        "Source": "No_AF2"}
        df.loc[len(df)] = motif_statistics
        df.to_csv(f"{cwd_path}motsasi_statistics.csv", index=False)
        return None
    print(f"{str(len(positives))} motifs in P2f")
    print(f"{str(len(negatives))} motifs in N2f")
    print(f"{str(len(remaining))} motifs in R2f")
    print(plotting_dict)

    # flexible iteration
    true_positives, positives, negatives, remaining, plotting_dict, motif_statistics, FinalMatrix = iteration_filter(motif_name, motif_label, motif_re, "flexible", screening_proteoma, true_positives, positives, negatives, remaining, Conservation_cutoff, SecStr_dict_cutoff, ExpBur_dict_cutoff, MinPos_cutoff, MaxPos_cutoff, top_GOTerms, clinvar_db_path, gnomad_db_path, uniref90_db_path, proteoma_list_path, sp_list_path, sp_fasta_path, trembl_fasta_path, tmp_path, pdbs_reparados_path, pdb_db_path, xml_db_path, sec_str_path, GO_db_path, sasa_thr_path, alphafold_human_path, scratch_path, cwd_path, plotting_dict, screen = screen)
    print(f"{str(len(positives))} motifs in P3f")
    print(f"{str(len(negatives))} motifs in N3f")
    print(f"{str(len(remaining))} motifs in R3f")
    print(plotting_dict)

    for positive in positives:
        print(positive.__dict__)

    iteration_process_plot(plotting_dict, motif_name, cwd_path, screen = screen)

    df.loc[len(df)] = motif_statistics

    df.to_csv(f"{cwd_path}motsasi_statistics.csv", index=False)

    variants_prediction(motif_name, true_positives, positives, FinalMatrix, cwd_path)

    negative_motifs(motif_name, negatives, cwd_path)

    return plotting_dict

#############################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Returns a list of match objects for the motif on a sequences database (default SwissProt)")
    parser.add_argument("Motif_name", action="store", help = "Name of the motif, e.g. DOC_MAPK_JIP1_4")
    parser.add_argument("Motif_re", action="store", help = "Motif regular expression, e.g. [RK]P[^P][^P]L.[LIVMF]")
    parser.add_argument("Motif_label", action="store", help = "Motif regular expression for plotting, e.g. RK.P.^P.^P.L.x.LIVMF")
    parser.add_argument("-cv_miss", "--clinvar_missense", action="store", help = "csv file with all clinvar missense variants", default="./ClinVar/ClinVar_missense_all_filtered.csv")
    parser.add_argument("-gad_miss", "--gnomad_missense", action="store", help = "csv file with all gnomad missense variants", default="./GnomAD/GnomAD_missense.csv")
    parser.add_argument("-elm_inst", "--elm_instances", action="store", help = "ELM Instances database where look for motifs", default="./ELM/elm_instances.tsv")
    parser.add_argument("-pdb_str", "--pdb_structures", action="store", help = "PDB database with all crystals", default="./PDB/zipped/")
    parser.add_argument("-tmp", "--tmp_path", action="store", help = "folder receiving foldx output", default="./tmp/")
    parser.add_argument("-rep_pdbs", "--repaired_pdbs", action="store", help = "folder with repaired PDBSs", default="./repaired_pdbs/")
    parser.add_argument("-uniref50", "--uniref50_clusters", action="store", help = "file with UniRef50 clusters for each UniProtID", default="./UniProt/UniRef50_sprot_h.csv")
    parser.add_argument("-sp_h_list", "--sp_h_list", action="store", help = "file with all UniProtIDs in SwissProt", default="./UniProt/uniprot_sprot_h.list")
    parser.add_argument("-sp", "--sp_fasta", action="store", help = "file with all fasta entries in SwissProt", default="./UniProt/uniprot_sprot.fasta.bgz")
    parser.add_argument("-sp_list", "--sp_list", action="store", help = "file with all UniProtIDs in SwissProt", default="./UniProt/uniprot_sprot.list")
    parser.add_argument("-trembl", "--trembl_fasta", action="store", help = "file with all fasta entries in TrEMBL", default="./UniProt/uniprot_trembl.fasta.bgz")
    parser.add_argument("-fasta", "--sp_h_fasta", action="store", help = "file with fasta entries for SwissProt", default="./UniProt/uniprot_sprot_h.fasta")
    parser.add_argument("-xml", "--sp_h_xml", action="store", help = "file with xml entries for SwissProt", default="./UniProt/uniprot_sprot_h.xml")
    parser.add_argument("-sec_est", "--secondary_structure", action="store", help = "folder with secondary structure predictions", default="./secondary_structures/")
    parser.add_argument("-exp_thr", "--exposure_thresholds", action="store", help = "exposure cut-offs", default="./SASA/Cortes_SASA.csv")
    parser.add_argument("-af_thr", "--af_threshold", action="store", help = "csv file with allele frequency for each UniProt, SwissProt entry gene", default="./GnomAD/uniprot_genes_freq_cutoffs.csv")
    parser.add_argument("-go", "--go_database", action="store", help = "file with codes for GO Terms", default="./GO/go.obo")
    parser.add_argument("-af_human", "--alphafold_human", action="store", help = "alphafold for the human proteome", default="./AlphaFold2_human_proteome/UP000005640_9606_HUMAN_v4/")
    parser.add_argument("-cwd", "--cwd", action="store", help = "current working directory path", default="./")
    parser.add_argument("-scratch", "--scratch_path", action="store", help = "scratch secondary structure predictor path", default="./scratch/SCRATCH-1D_1.3/bin/")
    parser.add_argument("-sp_omim_phen", "--sp_omim_phenotypes", action="store", help = "csv with phenotypes associated to UniProtIDs path", default="./OMIM/phenotype_uniprot.csv")

    args = parser.parse_args()

    if not os.path.exists(f"{args.cwd}Motifs/{args.Motif_name}"):
        os.system(f"mkdir {args.cwd}Motifs/{args.Motif_name}")

    if not os.path.exists(args.tmp_path):
        os.system(f"mkdir {args.tmp_path}")

    plotting_dict = iteration_process(args.Motif_name, args.Motif_re, args.Motif_label, args.clinvar_missense, args.gnomad_missense, args.af_threshold, args.uniref50_clusters, args.sp_h_list, args.sp_list, args.sp_fasta, args.trembl_fasta, args.tmp_path, args.sp_h_fasta, args.sp_h_xml, args.secondary_structure, args.go_database, args.elm_instances, args.pdb_structures, args.repaired_pdbs, args.alphafold_human, args.exposure_thresholds, args.cwd, args.scratch_path, args.sp_omim_phenotypes)
