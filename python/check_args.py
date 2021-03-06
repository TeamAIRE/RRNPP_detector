#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in gram-positive bacteria and bacteriophages genomes
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
    out_dir = os.path.join(parent_dir, 'rrnpp_detector_output')
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    return out_dir


def check_args(args):
    out_dir = check_out_dir(os.path.abspath(args.out_dir))
    # check user provided either a fna or a faa + an annotation file (either gff or feature table)
    if args.fna:
        if args.faa:
            sys.exit('error: Please choose between providing the fasta of the target genome (automated detection of ORFs) or of its protein sequences (pre-detected ORFs)')
        if args.gff:
            sys.exit('error: You cannot provide both the fasta of the target genome (automated detection of ORFs) and its annotations (pre-detected ORFs)')
        fasta, fasta_type = args.fna, 'fna'
    else:
        if not args.faa:
            sys.exit('error: Please provide a fasta of either the target genome (fna) or its predicted protein sequences (faa)')
        if args.gff and args.feature_tbl:
            sys.exit('error: Please choose between providing the genome annotation in gff format or in NCBI assembly feature table format')
        elif args.gff:
            check_file_exists(args.gff)
            check_file_extension(args.gff, 'gff')
        elif args.feature_tbl:
            check_file_exists(args.feature_tbl)
        else:
            sys.exit('error: If you provide the fasta of the protein sequences of the target genome, you also have to provide its annotations (in gff or NCBI assembly feature table formats)')
        fasta, fasta_type = args.faa, 'faa'
    check_fasta(fasta, fasta_type)
    try:
        int(args.cpu)
    except:
        sys.exit('error: Please provide an integer for the number of cpu to be used')
    return out_dir, fasta, fasta_type
            
 