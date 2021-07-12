#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in gram-positive bacteria and bacteriophages genomes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import os
import pandas
import subprocess


# apply thresholds and keep only the best homolog for each protein   
def filter_blastp_results(blastp_results, blastp_min_pident, blastp_min_pcover, blastp_max_evalue):
    df = pandas.read_csv(blastp_results, sep='\t', header=None)
    df.columns = ['homologous_ref_receptor', 'target_name', 'qlen', 'slen', 'perc_identity', 'blastp_evalue', 'bitscore', 'qstart', 'qend', 'sstart', 'send']
    query_coverage = round((df['qend'] - df['qstart'] + 1) / df['qlen'] * 100, 1)
    target_coverage = round((df['send'] - df['sstart'] + 1) / df['slen'] * 100, 1)
    df['perc_cover'] = pandas.concat([query_coverage, target_coverage], axis=1).apply(min, axis=1)
    filtered_df = df.loc[(df['perc_identity'] >= blastp_min_pident) & (df['perc_cover'] >= blastp_min_pcover) & (df['blastp_evalue'] <= blastp_max_evalue)]
    filtered_df = filtered_df.sort_values(by = 'blastp_evalue', ascending = True, kind = 'mergesort').drop_duplicates(subset=['target_name'])
    return filtered_df[['target_name', 'homologous_ref_receptor', 'blastp_evalue', 'perc_identity', 'perc_cover']]


def blastp(target_faa, out_dir, cpu, rrnpp_detector_dir, blastp_min_pident, blastp_min_pcover, blastp_max_evalue):
    # makeblastdb
    database = os.path.join(out_dir, 'candidate_receptors_blastdb')
    makeblastdb_args = ['makeblastdb', '-dbtype', 'prot', '-in', target_faa, '-out', database]
    subprocess.run(makeblastdb_args, check = True)
    
    # blastp
    query_faa = os.path.join(rrnpp_detector_dir, 'reference_receptors', 'reference_receptors.faa')
    blastp_results = os.path.join(out_dir, 'reference_receptors_vs_candidate_receptors_blastp_out.tsv')
    blastp_args = ['blastp', '-query', query_faa, '-db', database, '-out', blastp_results, '-num_threads', cpu, 
                   '-max_target_seqs', '1000000', '-outfmt', '6 qseqid sseqid qlen slen pident evalue bitscore qstart qend sstart send']
    subprocess.run(blastp_args, check = True)

    # clean blast
    df = filter_blastp_results(blastp_results, blastp_min_pident, blastp_min_pcover, blastp_max_evalue)
    return df
    
            