#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2025

@author: goldenpolaco
"""
import math, sys, getopt
from Bio import ExPASy
from Bio import SeqIO
from lxml import etree
import requests
from io import StringIO, BytesIO
from tqdm import tqdm
import os
import pandas as pd
import sys
import urllib
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# numarray imported below
PSEUDOCOUNT = .0000001

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"] 

# dictionary to map from amino acid to its row/column in a similarity matrix
aa_to_index = {}
for i, aa in enumerate(amino_acids):
    aa_to_index[aa] = i

################################################################################
# Frequency Count and Gap Penalty
################################################################################

def weighted_freq_count_pseudocount(col, seq_weights, pc_amount):
    """ Return the weighted frequency count for a column--with pseudocount."""

    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1. for i in range(len(col))]

    aa_num = 0
    freq_counts = [pc_amount for i in range(len(amino_acids))] # in order defined by amino_acids

    for aa in amino_acids:
        for j in range(len(col)):
            if col[j] == aa:
                freq_counts[aa_num] += 1 * seq_weights[j]

        aa_num += 1

    for j in range(len(freq_counts)):
        freq_counts[j] = freq_counts[j] / (sum(seq_weights) + len(amino_acids) * pc_amount)

    return freq_counts


def weighted_gap_penalty(col, seq_weights):
    """ Calculate the simple gap penalty multiplier for the column. If the 
    sequences are weighted, the gaps, when penalized, are weighted 
    accordingly. """

    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.] * len(col)
    
    gap_sum = 0.
    for i in range(len(col)):
        if col[i] == '-':
            gap_sum += seq_weights[i]

    return 1 - (gap_sum / sum(seq_weights))


def gap_percentage(col):
    """Return the percentage of gaps in col."""
    num_gaps = 0.

    for aa in col:
        if aa == '-': num_gaps += 1

    return num_gaps / len(col)

################################################################################
# Jensen-Shannon Divergence
################################################################################

def js_divergence(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
    """
    Return the Jensen-Shannon Divergence for the column with the background
    distribution bg_distr. sim_matrix is ignored. JSD is the default method
    """

    distr = bg_distr[:]

    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    # if background distrubtion lacks a gap count, remove fc gap count
    if len(distr) == 20: 
        new_fc = fc[:-1]
        s = sum(new_fc)
        for i in range(len(new_fc)):
            new_fc[i] = new_fc[i] / s
        fc = new_fc

    if len(fc) != len(distr): return -1

    # make r distriubtion
    r = []
    for i in range(len(fc)):
        r.append(.5 * fc[i] + .5 * distr[i])

    d = 0.
    for i in range(len(fc)):
        if r[i] != 0.0:
            if fc[i] == 0.0:
                d += distr[i] * math.log(distr[i]/r[i], 2)
            elif distr[i] == 0.0:
                d += fc[i] * math.log(fc[i]/r[i], 2) 
            else:
                d += fc[i] * math.log(fc[i]/r[i], 2) + distr[i] * math.log(distr[i]/r[i], 2)

    # d /= 2 * math.log(len(fc))
    d /= 2

    if gap_penalty == 1: 
        return d * weighted_gap_penalty(col, seq_weights)
    else: 
        return d



################################################################################
# Window Score
################################################################################

def window_score(scores, window_len, lam=.5):
    """
    This function takes a list of scores and a length and transforms them 
    so that each position is a weighted average of the surrounding positions. 
    Positions with scores less than zero are not changed and are ignored in the 
    calculation. Here window_len is interpreted to mean window_len residues on 
    either side of the current residue
    """

    w_scores = scores[:]

    for i in range(window_len, len(scores) - window_len):
        if scores[i] < 0: 
            continue

        sum = 0.
        num_terms = 0.
        for j in range(i - window_len, i + window_len + 1):
            if i != j and scores[j] >= 0:
                num_terms += 1
                sum += scores[j]

        if num_terms > 0:
            w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[i]

    return w_scores


def calc_z_scores(scores, score_cutoff):
    """
    Calculates the z-scores for a set of scores. Scores below
    score_cutoff are not included
    """

    average = 0.
    std_dev = 0.
    z_scores = []
    num_scores = 0

    for s in scores:
        if s > score_cutoff:
            average += s
            num_scores += 1
    if num_scores != 0:
        average /= num_scores

    for s in scores:
        if s > score_cutoff:
            std_dev += ((s - average)**2) / num_scores
    std_dev = math.sqrt(std_dev)

    for s in scores:
        if s > score_cutoff and std_dev != 0:
            z_scores.append((s-average)/std_dev)
        else:
            z_scores.append(-1000.0)

    return z_scores


################################################################################
################################################################################
################################################################################
#  END CONSERVATION SCORES
################################################################################
################################################################################
################################################################################

def read_scoring_matrix(sm_file):
    """
    Read in a scoring matrix from a file, e.g., blosum80.bla, and return it
    as an array
    """
    aa_index = 0
    first_line = 1
    row = []
    list_sm = [] # hold the matrix in list form

    try:
        matrix_file = open(sm_file, 'r')

        for line in matrix_file:

            if line[0] != '#' and first_line:
                first_line = 0
                if len(amino_acids) == 0:
                    for c in line.split():
                        aa_to_index[string.lower(c)] = aa_index
                        amino_acids.append(string.lower(c))
                        aa_index += 1

            elif line[0] != '#' and first_line == 0:
                if len(line) > 1:
                    row = line.split()
                    list_sm.append(row)

    except IOError:
        print("Could not load similarity matrix: %s. Using identity matrix..." % sm_file)
        return identity(20)
        
    # if matrix is stored in lower tri form, copy to upper
    if len(list_sm[0]) < 20:
        for i in range(0,19):
            for j in range(i+1, 20):
                list_sm[i].append(list_sm[j][i])

    for i in range(len(list_sm)):
        for j in range(len(list_sm[i])):
            list_sm[i][j] = float(list_sm[i][j])

    return list_sm
    #sim_matrix = array(list_sm, type=Float32)
    #return sim_matrix

def calculate_sequence_weights(msa):
    """
    Calculate the sequence weights using the Henikoff '94 method
    for the given msa
    """

    seq_weights = [0. for a in range(len(msa))]
    
    for i in range(len(msa[0])):

        freq_counts = [0 for b in range(len(amino_acids))]
        col = []
        
        for j in range(len(msa)):
            if msa[j][i] != '-': # ignore gaps
                freq_counts[aa_to_index[msa[j][i]]] += 1

        num_observed_types = 0
        for j in range(len(freq_counts)):
            if freq_counts[j] > 0:
                num_observed_types +=1

        for j in range(len(msa)):
            d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
            if d > 0:
                seq_weights[j] += 1. / d

    for w in range(len(seq_weights)):
            seq_weights[w] /= len(msa[0])

    return seq_weights


def load_sequence_weights(fname):
    """
    Read in a sequence weight file f and create sequence weight list. 
    The weights are in the same order as the sequences each on a new line
    """
    
    seq_weights = []

    try:
        f = open(fname)

        for line in f:
            l = line.split()
            if line[0] != '#' and len(l) == 2:
                seq_weights.append(float(l[1]))

    except IOError:
        #print ("No sequence weights. Can't find: " + fname)
        pass

    return seq_weights

def get_column(col_num, alignment):
    """
    Return the col_num column of alignment as a list.
    """
    col = []
    for seq in alignment:
        if col_num < len(seq):
            col.append(seq[col_num])
    return col

def get_distribution_from_file(fname):
    """
    Read an amino acid distribution from a file. The probabilities should
    be on a single line separated by whitespace in alphabetical order as in 
    amino_acids above. # is the comment character
    """

    distribution = []
    try:
        f = open(fname)
        for line in f:
            if line[0] == '#':
                continue
            line = line[:-1]
            distribution = line.split()
            distribution = map(float, distribution)

            
    except IOError:
        print("Using default (BLOSUM62) background.")
        return []

    # use a range to be flexible about round off
    if .997 > sum(distribution) or sum(distribution) > 1.003:
        print("Distribution does not sum to 1. Using default (BLOSUM62) background.")
        print(sum(distribution))
        return []

    return distribution


def read_fasta_alignment(filename):
    """
    Read in the alignment stored in the FASTA file, filename. Return two
    lists: the identifiers and sequences
    """

    f = open(filename)

    names = []
    alignment = []
    cur_seq = ''

    for line in f:
        line = line[:-1]
        if len(line) == 0:
            continue
        if line[0] == ';':
            continue
        if line[0] == '>':
            names.append(line[1:].replace('\r', ''))

            if cur_seq != '':
                cur_seq = cur_seq.upper()
                for i, aa in enumerate(cur_seq):
                    if aa not in iupac_alphabet:
                        cur_seq = cur_seq.replace(aa, '-')
                alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
                cur_seq = ''
        elif line[0] in iupac_alphabet:
            cur_seq += line.replace('\r', '')

    # add the last sequence
    cur_seq = cur_seq.upper()
    for i, aa in enumerate(cur_seq):
        if aa not in iupac_alphabet:
            cur_seq = cur_seq.replace(aa, '-')
    alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))

    return names, alignment

def read_clustal_alignment(filename):
    """
    Read in the alignment stored in the CLUSTAL file, filename. Return
    two lists: the names and sequences
    """

    names = []
    alignment = []

    f = open(filename)

    for line in f:
        line = line[:-1]
        if len(line) == 0:
            continue
        if '*' in line:
            continue

        if 'CLUSTAL' in line:
            continue

        t = line.split()

        if len(t) == 2 and t[1][0] in iupac_alphabet:
            if t[0] not in names:
                names.append(t[0])
                alignment.append(t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X', '-').replace('\r', ''))
            else:
                alignment[names.index(t[0])] += t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X','-').replace('\r', '')
                   
    return names, alignment

################################################################################
# Begin execution
################################################################################

#UniProtID = str(input("Ingrese el UniProtID de la proteína consultada: "))
UniProtID = sys.argv[1]
tmp_path = sys.argv[2]

# BLOSUM62 background distribution
blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]

# set defaults
window_size = 3 # 0 = no window
win_lam = .5 # for window method linear combination
outfile_name = ""
s_matrix_file = "./blosum_matrix/blosum62.bla"
bg_distribution = blosum_background_distr[:]
scoring_function = js_divergence
use_seq_weights = True
background_name = 'blosum62'
gap_cutoff = .3
use_gap_penalty = 1
seq_specific_output = UniProtID # name of sequence if True
normalize_scores = False

align_file = f"{tmp_path}mafft_{UniProtID}.fasta"
align_suffix = align_file.split('.')[-1]

s_matrix = read_scoring_matrix(s_matrix_file)

names = []
alignment = []
seq_weights = []

try:
    names, alignment = read_clustal_alignment(align_file)
    if names == []:
        names, alignment = read_fasta_alignment(align_file)
except IOError:
    print("Could not find %s. Exiting..." % align_file)
    sys.exit(1)


if len(alignment) != len(names) or alignment == []:
    print("Unable to parse alignment.\n")
    sys.exit(1)

seq_len = len(alignment[0])
for i, seq in enumerate(alignment):
    if len(seq) != seq_len:
        print("ERROR: Sequences of different lengths: %s (%d) != %s (%d).\n" % (names[0], seq_len, names[i], len(seq)))
        sys.exit(1)


if use_seq_weights:
    seq_weights = load_sequence_weights(align_file.replace('.%s' % align_suffix, '.weights'))
    if seq_weights == []:
        seq_weights = calculate_sequence_weights(alignment)

if len(seq_weights) != len(alignment): seq_weights = [1. for i in range(len(alignment))]

# handle print of output relative to specific sequence
ref_seq_num = None
if seq_specific_output and seq_specific_output not in names:
    print("Sequence %s not found in alignment. Using default output format...\n" % seq_specific_output)
    seq_specific_output = 0
elif seq_specific_output in names:
    ref_seq_num = names.index(seq_specific_output)

# calculate scores
scores = []
for i in range(len(alignment[0])):
    col = get_column(i, alignment)

    if len(col) == len(alignment):
        if gap_percentage(col) <= gap_cutoff:
            scores.append(scoring_function(col, s_matrix, bg_distribution, seq_weights, use_gap_penalty))
        else:
            scores.append(-1000.)

if window_size > 0:
    scores = window_score(scores, window_size, win_lam)

if normalize_scores:
    scores = calc_z_scores(scores, -999)

############## agregado por Franco ########################

scores_interes = []
aas = []
count = 1
aas_number = []

if seq_specific_output:
    for i in range(len(alignment[0])):
        cur_aa = get_column(i, alignment)[ref_seq_num]
        if cur_aa == '-':
            continue
        else:
            scores_interes.append(scores[i])
            aas.append(cur_aa)
            aas_number.append(count)
            count+=1
            
# print(len(scores_interes))
# print(aas)
# print(aas_number)

b = {"Residue_number": aas_number, "Residue": aas, "Bitscores": scores_interes}

bb = pd.DataFrame(b, columns = ["Residue_number", "Residue", "Bitscores"])

#print(bb)
bb.to_csv(f"{tmp_path}Conservacion_jsdiv_{UniProtID}.tsv", sep = "\t", index=False)
