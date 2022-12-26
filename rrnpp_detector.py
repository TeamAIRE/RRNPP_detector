#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-type quorum sensing systems in chromosomes, plasmids and phages of Firmicutes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import argparse
import os
import sys
import rrnpp_detector
from rrnpp_detector.check_args import check_args
from rrnpp_detector.preprocessing import load_parameters, create_directories
from rrnpp_detector.integrate_annotations import get_annotation_config
from rrnpp_detector.launcher import launcher
    

def main():
    # Read arguments
    parser = argparse.ArgumentParser(description='RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems '
                                     'in chromosomes, plasmids and bacteriophages of Firmicutes')
    parser.add_argument('--version', action='version', help='print version number and exit.', version=rrnpp_detector.__version__)
    # Input files
    parser.add_argument('-o', dest='out_dir', type=str, default='./', help='path to output directory (default is current directory)')
    parser.add_argument('--fna', dest='fna', type=str, help='path to the fasta of the target genome(s) (will run Prodigal to detect CDSs if faa not provided)')
    parser.add_argument('--faa', dest='faa', type=str, help='path to fasta of the protein sequences of the target genome(s) '
                        '(requires additional --gff or --ft option)')
    parser.add_argument('--gff', dest='gff', type=str, help='path to the annotations of the target genome(s) in gff')
    parser.add_argument('--ft', dest='feature_tbl', type=str, help='path to the annotations of the target genome(s) in the NCBI_assembly feature_table format')
    # Run options
    parser.add_argument('--cpu', dest='cpu', type=str, default='1', help='number of cpu to use (default is 1)')
    parser.add_argument('--chunk_size', dest='chunk_size', help='nb target genomes to be processed altogether to preserve RAM usage '
                        '(e.g. if --chunk_size 100, then the program will process 100 genomes by 100 genomes instead of all together)')
    parser.add_argument('--keep_working_dir', action='store_true', help='keep the directory of intermediate files')
    # Search options
    parser.add_argument('--min_pl',  dest='min_propeptide_len', default='10',  help='minimal propeptide length (default=10)')
    parser.add_argument('--max_pl',  dest='max_propeptide_len', default='100',  help='maximal propeptide length (default=100)')
    parser.add_argument('--min_rl',  dest='min_receptor_len',   default='250', help='minimal receptor length (default=250)')
    parser.add_argument('--max_rl',  dest='max_receptor_len',   default='500', help='maximal receptor length (default=500)')
    parser.add_argument('--min_igd', dest='min_intergen_dist',  default='-60', help='minimal intergenic distance (default=-60)')
    parser.add_argument('--max_igd', dest='max_intergen_dist',  default='400', help='maximal intergenic distance (default=400)')
    parser.add_argument('--start_codons', dest='start_codons', default='ATG', help='comma-separated list of start codons to consider for ORF calling (default=ATG)')
    parser.add_argument('--rbs_bins', dest='rbs_bins', default='27,24,23,22,20,19,16,15,14,13,12,6', help='comma-separated list of Prodigal\'s RBS bins to consider for ORF calling '
                        '(default=27,24,23,22,20,19,16,15,14,13,12,6), '
                        'to by bypass the filter, use --rbs_bins 27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0')
    parser.add_argument('--expand_to_homologs', action='store_true', help='use detected systems as seeds to detect putative homologous systems missed by RRNPP_detector')
    parser.add_argument('--tprpred', dest='tprpred', action='store_true', help='run tprpred in addition to hmmsearch for TPR motifs detection')
    parser.add_argument('--predisi', dest='predisi', action='store_true', help='run predisi in addition to signalp for the detection of a secretion-tag within propeptides')

    args = parser.parse_args()
    rrnpp_detector_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    current_dir = os.path.abspath(os.getcwd())      
 
    # Check args and files provided by the user are OK
    out_dir, faa, fna = check_args(args)
    print('Running rrnpp_detector v' + rrnpp_detector.__version__)
    
    # Parameters
    parameters = load_parameters(args)

    # Create output, working, fasta directories
    out_dir, working_dir, fasta_dir = create_directories(out_dir)
    
    # Get annotation config
    faa, annotation_file, annotation_format = get_annotation_config(faa, fna, out_dir, args)
        
    # Run rrnpp_detector
    launcher(args, parameters, rrnpp_detector_dir, current_dir, out_dir, working_dir, fasta_dir, 
             fna, faa, annotation_file, annotation_format)
      
if __name__=='__main__': 
    main()

