#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-type quorum sensing systems in chromosomes, plasmids and phages of Firmicutes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import os
import re
import sys
import shutil


def check_file_exists(file):
    if not os.path.exists(file):
        sys.exit('error: Unable to locate \'%s\'' % file)

  
def check_file_extension(file, file_type):
    if file_type == 'fna':
        if not re.match(r'.*?\.f(n)?a(sta)?$', file):
            sys.exit('error: File \'%s\' has an incorrect extension: \'.fa\', \'.fna\' or \'.fasta\' are expected for the fasta of the target genome.' % file)
    elif file_type == 'faa':
        if not re.match(r'.*?\.f(a)?a(sta)?$', file):
            sys.exit('error: File \'%s\' has an incorrect extension: \'.fa\', \'.faa\' or \'.fasta\' are expected for the fasta of the protein sequences of the target genome.' % file)
    elif file_type == 'gff':
        if not re.match(r'.*?\.gff(3)?', file):
            sys.exit('error: File \'%s\' has an incorrect extension: \'.gff\', \'.gff3\' are expected for the genome annotation.' % file)
    
              
def check_fasta(file, fasta_type):
    check_file_exists(file)
    check_file_extension(file, fasta_type)
    with open(file, mode = 'r') as f:  
        line1 = f.readline()
        # check first character of the fasta is '>'
        if not line1.startswith('>'):
            sys.exit('error: Fasta \'%s\' does not begin with the expected \'>\' character' % file)
        # check second line of the fasta correspond either to either a nucleotide sequence (fna) or a protein sequence (faa)
        line2 = f.readline().strip('\n')
        if fasta_type == 'fna':
            if not re.match(r'^[ACTGX -*]+$', line2):
                sys.exit('error: Second line of the fasta \'%s\' does not correspond to a nucleotide sequence' % file)
        elif fasta_type == 'faa':
            if not re.match(r'^[A-Z -*]+$', line2):
                sys.exit('error: Second line of the fasta \'%s\' does not correspond to a protein sequence' % file)
            elif re.match(r'^[ACTGX -*]+$', line2):
                sys.exit('error: Second line of the fasta \'%s\' seems to correspond to a nucleotide sequence and not to a protein sequence' % file)
    f.close()


def check_out_dir(parent_dir):
    # create a directory 'rrnpp_detector_output' in the output directory provided by the user
    if not os.path.exists(parent_dir):
        sys.exit('error: Unable to locate \'%s\'. Please provide a correct path to the output directory' % parent_dir)
        

def rbs_bins_to_list(rbs_bins_arg):
    try:
        bin_list = [int(rbs_bin) for rbs_bin in rbs_bins_arg.split(',')]
    except:
        sys.exit('error: argument --rbs_bins only accepts a comma-separated list of integers (from 0 to 27)')
    if min(bin_list) < 0:
        sys.exit('error: minimal RBS bin allowed is 0 (no RBS)')
    if max(bin_list) > 27:
        sys.exit('error: maximial RBS bin allowed is 27')
    return bin_list


def start_codons_to_list(start_codons_arg):
    try:
        start_codons_list = [start_codon.upper() for start_codon in start_codons_arg.split(',')]
    except:
        sys.exit('error: argument --start_codons only accepts a comma-separated list of start codons')
    for codon in start_codons_list:
        if not re.match(r'^[A-Z]{3}$', codon):
            sys.exit('error: %s is not a codon' % codon)
    return start_codons_list


def check_args(args):
    # combinations allowed:
    # - fna                     (Prodigal + Small orf search)
    # - fna + faa + annotations (Small orf search)
    # - faa + annotations       (No small orf search)
    out_dir = os.path.abspath(args.out_dir)
    check_out_dir(out_dir)
    faa = None
    fna = None
    if args.faa:
        if args.gff and args.feature_tbl:
            sys.exit('error: Please choose between providing the genome annotation in gff format or in NCBI assembly feature table format')
        elif args.gff:
            check_file_exists(args.gff)
            check_file_extension(args.gff, 'gff')
        elif args.feature_tbl:
            check_file_exists(args.feature_tbl)        
        else:
            sys.exit('error: If you provide the fasta of the protein sequences of the target genome, '
                     'you also have to provide its annotations (in gff or NCBI assembly feature table formats)')
        faa = args.faa
        check_fasta(faa, 'faa')
        if args.fna:
            fna = args.fna
            check_fasta(fna, 'fna')
    elif args.fna:
        if args.gff or args.feature_tbl:
            sys.exit('error: You cannot provide the annotations of the genome(s) without the fasta of the protein sequences')
        fna = args.fna
        check_fasta(fna, 'fna')
    else:
        sys.exit('error: Please provide a fasta of either the target genome (fna) or its predicted protein sequences (faa)')
    try:
        int(args.cpu)
    except:
        sys.exit('error: Please provide an integer for the number of cpu to be used')
    args.rbs_bins = rbs_bins_to_list(args.rbs_bins)
    args.start_codons = start_codons_to_list(args.start_codons)
    return out_dir, faa, fna
            
 