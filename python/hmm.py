#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in gram-positive bacteria and bacteriophages genomes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import pandas
import re
import subprocess


# because the tblout output of hmmsearch is not a tsv but a messy ' '-separated file
def tblout_to_tsv(tblout, tsv):
    with open(tsv, mode='w') as outfile:
        header = 'target_name\tt_accession\thmm_name\thmm_accession\tE-value\tscore\tbias\tE-value_best_domain\tscore_best_domain\tbias_best_domain\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tdescription_of_target\n'
        outfile.write(header)    
        with open(tblout, mode='r') as infile:
            for line in infile:
                if not line.startswith('#'):
                    newline = re.sub('  *', '\t', line, count = 18)
                    outfile.write(newline)
        infile.close()
    outfile.close()


# for each protein matched (evalue <= max_evalue), keep only the info of the hmm that gave rise to the best evalue
def filter_hmmsearch_results(tsv, max_evalue):
    df = pandas.read_csv(tsv, sep='\t')
    df = df.loc[df['E-value'] <= max_evalue]
    df = df.sort_values(by='E-value', ascending=True, kind='mergesort').drop_duplicates(subset=['target_name']) 
    return df     


def hmmsearch(faa, hmm_library, stdout, tblout, cpu, max_evalue):
    hmmsearch_args = ['hmmsearch', '--cpu', cpu, '--noali', '--tblout', tblout,
                      '-o', stdout, hmm_library, faa]
    subprocess.run(hmmsearch_args, check = True)
    
    tsv = tblout[:-4] + '.tsv'
    tblout_to_tsv(tblout, tsv)
    results_df = filter_hmmsearch_results(tsv, max_evalue)
    retained_receptors = results_df['target_name'].tolist()

    return retained_receptors, results_df[['target_name', 'hmm_name', 'hmm_accession', 'E-value']]  
    
            