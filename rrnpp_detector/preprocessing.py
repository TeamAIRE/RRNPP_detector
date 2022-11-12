#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-type quorum sensing systems in chromosomes, plasmids and phages of Firmicutes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import os
import re


def load_parameters(args):
    parameters = dict()
    parameters['propeptide_len_boundaries'] = [int(args.min_propeptide_len), int(args.max_propeptide_len)]
    parameters['receptor_len_boundaries'] = [int(args.min_receptor_len), int(args.max_receptor_len)]
    parameters['intergenic_distance_boundaries'] = [int(args.min_intergen_dist), int(args.max_intergen_dist)]
    parameters['hmmsearch_max_evalue'] = 1E-4
    parameters['tprpred_min_probab'] = 10
    parameters['predisi_min_likelihood'] = 0.3
    parameters['blastp_max_evalue'] = 1E-5
    parameters['blastp_min_pident'] = 20
    parameters['blastp_min_pcover'] = 60
    parameters['blastp_max_evalue_iterative_search_receptors'] = 1E-5
    parameters['blastp_min_pident_iterative_search_receptors'] = 25
    parameters['blastp_min_pcover_iterative_search_receptors'] = 70
    parameters['blastp_max_evalue_iterative_search_propeptides'] = 0.05
    parameters['blastp_min_pident_iterative_search_propeptides'] = 25
    parameters['blastp_min_pcover_iterative_search_propeptides'] = 50
    # parameters['orfipy_start_codons'] = ['ATG', 'GTG', 'TTG']
    parameters['orfipy_start_codons'] = ['ATG']
    # parameters['allowed_rbs_bins'] = [i for i in range(27, -1, -1)]
    parameters['allowed_rbs_bins'] = [27, 24, 23, 22, 20, 19, 16, 15, 14, 13, 12, 6]
    parameters['most_used_rbs_bins'] = [27, 24, 22, 19, 16, 15, 14, 13, 12, 6] # defined according to Omotajo et al.
    parameters['strict_contexts'] = ['(====> =>)', '(<= <====)', '(<= ====>)', '(<==== =>)'] # divergent + co-directional with peptide downstream of receptor
    parameters['strict_predictors'] = ['SP(Sec/SPI)', 'SHP_propeptide']
    parameters['relaxed_predictors'] = ['SP(Sec/SPI)', 'SHP_propeptide', 'LIPO(Sec/SPII)', 'TAT(Tat/SPI)', 'PrediSi_pos', 'PrediSi_neg_(but_high_score)']
    return parameters


def create_directories(parent_dir):
    # Output directory
    out_dir = prefix = os.path.join(parent_dir, 'rrnpp_detector_output')
    k = 0
    while os.path.exists(out_dir):
        k += 1
        out_dir = prefix + '_' + str(k)
    os.makedirs(out_dir)
    # Working directory
    working_dir = os.path.join(out_dir, 'working_dir')
    os.makedirs(working_dir)
    # Fasta directory
    fasta_dir = os.path.join(out_dir, 'protein_sequences')
    os.makedirs(fasta_dir)
    return out_dir, working_dir, fasta_dir
    

def load_rbs_regex():
    rbs_dict = dict()
    for i in range(0, 30):
        rbs_dict[i] = dict()
        
    rbs_dict[0]["bin"] = 27
    rbs_dict[0]["length"] = 6
    rbs_dict[0]["regex"] = re.compile(r"AGGAGG.{5,10}$")
    
    rbs_dict[1]["bin"] = 26
    rbs_dict[1]["length"] = 6
    rbs_dict[1]["regex"] = re.compile(r"AGGAGG.{3,4}$")

    rbs_dict[2]["bin"] = 25
    rbs_dict[2]["length"] = 6
    rbs_dict[2]["regex"] = re.compile(r"AGGAGG.{11,12}$")

    rbs_dict[3]["bin"] = 24
    rbs_dict[3]["length"] = 5
    rbs_dict[3]["regex"] = re.compile(r"GGAGG.{5,10}$")

    rbs_dict[4]["bin"] = 23
    rbs_dict[4]["length"] = 5
    rbs_dict[4]["regex"] = re.compile(r"GGAGG.{3,4}$")

    rbs_dict[5]["bin"] = 22
    rbs_dict[5]["length"] = 5
    rbs_dict[5]["regex"] = re.compile(r"AGGAG.{5,10}$")

    rbs_dict[6]["bin"] = 21
    rbs_dict[6]["length"] = 5
    rbs_dict[6]["regex"] = re.compile(r"AGGAG.{3,4}$")

    rbs_dict[7]["bin"] = 20
    rbs_dict[7]["length"] = 5
    rbs_dict[7]["regex"] = re.compile(r"(AGGAG|GGAGG).{11,12}$")

    rbs_dict[8]["bin"] = 19
    rbs_dict[8]["length"] = 6
    rbs_dict[8]["regex"] = re.compile(r"(AG.{1}AGG|AGG.{1}GG).{5,10}$")

    rbs_dict[9]["bin"] = 18
    rbs_dict[9]["length"] = 6
    rbs_dict[9]["regex"] = re.compile(r"(AG.{1}AGG|AGG.{1}GG).{3,4}$")

    rbs_dict[10]["bin"] = 17
    rbs_dict[10]["length"] = 6
    rbs_dict[10]["regex"] = re.compile(r"(AG.{1}AGG|AGG.{1}GG).{11,12}$")

    rbs_dict[11]["bin"] = 16
    rbs_dict[11]["length"] = 4
    rbs_dict[11]["regex"] = re.compile(r"(GGAG|GAGG).{5,10}$")

    rbs_dict[12]["bin"] = 15
    rbs_dict[12]["length"] = 4
    rbs_dict[12]["regex"] = re.compile(r"AGGA.{5,10}$")

    rbs_dict[13]["bin"] = 14
    rbs_dict[13]["length"] = 5
    rbs_dict[13]["regex"] = re.compile(r"GG.{1}GG.{5,10}$")

    rbs_dict[14]["bin"] = 13
    rbs_dict[14]["length"] = 3
    rbs_dict[14]["regex"] = re.compile(r"(GGA|GAG|AGG).{5,10}$")

    rbs_dict[15]["bin"] = 12
    rbs_dict[15]["length"] = 4
    rbs_dict[15]["regex"] = re.compile(r"(AGGA|GGAG|GAGG).{11,12}$")

    rbs_dict[16]["bin"] = 11
    rbs_dict[16]["length"] = 4
    rbs_dict[16]["regex"] = re.compile(r"(AGGA|GGAG|GAGG).{3,4}$")

    rbs_dict[17]["bin"] = 10
    rbs_dict[17]["length"] = 6
    rbs_dict[17]["regex"] = re.compile(r"AGGAGG.{13,15}$")

    rbs_dict[18]["bin"] = 10
    rbs_dict[18]["length"] = 5
    rbs_dict[18]["regex"] = re.compile(r"(AGGAG|GGAGG).{13,15}$")

    rbs_dict[19]["bin"] = 9
    rbs_dict[19]["length"] = 5
    rbs_dict[19]["regex"] = re.compile(r"AG.{1}AG.{5,10}$")

    rbs_dict[20]["bin"] = 8
    rbs_dict[20]["length"] = 5
    rbs_dict[20]["regex"] = re.compile(r"GG.{1}GG.{3,4}$")

    rbs_dict[21]["bin"] = 7
    rbs_dict[21]["length"] = 5
    rbs_dict[21]["regex"] = re.compile(r"GG.{1}GG.{11,12}$")

    rbs_dict[22]["bin"] = 6
    rbs_dict[22]["length"] = 3
    rbs_dict[22]["regex"] = re.compile(r"(GGA|GAG|AGG).{11,12}$")

    rbs_dict[23]["bin"] = 5
    rbs_dict[23]["length"] = 5
    rbs_dict[23]["regex"] = re.compile(r"AG.{1}AG.{3,4}$")

    rbs_dict[24]["bin"] = 4
    rbs_dict[24]["length"] = 5
    rbs_dict[24]["regex"] = re.compile(r"AG.{1}AG.{11,12}$")

    rbs_dict[25]["bin"] = 3
    rbs_dict[25]["length"] = 6
    rbs_dict[25]["regex"] = re.compile(r"(AG.{1}AGG|AGG.{1}GG).{13,15}$")

    rbs_dict[26]["bin"] = 3
    rbs_dict[26]["length"] = 4
    rbs_dict[26]["regex"] = re.compile(r"(AGGA|GGAG|GAGG).{13,15}$")

    rbs_dict[27]["bin"] = 2
    rbs_dict[27]["length"] = 5
    rbs_dict[27]["regex"] = re.compile(r"(AG.{1}AG|GG.{1}GG).{13,15}$")

    rbs_dict[28]["bin"] = 2
    rbs_dict[28]["length"] = 3
    rbs_dict[28]["regex"] = re.compile(r"(GGA|GAG|AGG).{13,15}$")

    rbs_dict[29]["bin"] = 1
    rbs_dict[29]["length"] = 3
    rbs_dict[29]["regex"] = re.compile(r"(GGA|GAG|AGG).{3,4}$")
    return rbs_dict


def load_translation_table():
    # "name": "Bacterial, Archaeal and Plant Plastid (transl_table=11)",
    # "start": ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
    # "stop": ["TAA","TAG","TGA"]
    translation_table={
        "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", 
        "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
        "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
        "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
        "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
        "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
        "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
        "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
        "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
        "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
        "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
        "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}
    return translation_table
