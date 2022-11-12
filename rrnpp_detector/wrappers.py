#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-type quorum sensing systems in chromosomes, plasmids and phages of Firmicutes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import os
import pandas
import re
import shutil
import subprocess


#######################################
# Prodigal
#######################################
def prodigal(fna, out_dir):
    # produce a faa and a gff from a fna
    gff = os.path.join(out_dir, 'predicted_ORFs.gff')
    faa = os.path.join(out_dir, 'predicted_proteins.faa')
    prodigal_args = ['prodigal', '-i', fna, '-f', 'gff', '-o', gff, '-a', faa]
    subprocess.run(prodigal_args, check = True)
    print('  Prodigal done!')
    return gff, faa


#######################################
# hmmsearch
#######################################
# because the tblout output of hmmsearch is not a tsv but a messy ' '-separated file
def tblout_to_tsv(tblout, tsv, header, nb_field_separators):
    with open(tsv, mode='w') as outfile:
        if header:
            outfile.write(header)     
        with open(tblout, mode='r') as infile:
            for line in infile:
                if not line.startswith('#'):
                    newline = re.sub('  *', '\t', line, count = nb_field_separators)
                    outfile.write(newline)
        infile.close()
    outfile.close()


# for each protein matched (evalue <= max_evalue), keep only the info of the hmm that gave rise to the best evalue
def filter_hmmsearch_results(tsv, max_evalue):
    useful_column_names = ['target_name', 'hmm_name', 'hmm_accession', 'E-value', 'score']
    column_types = {'target_name': str, 'hmm_name': str, 'hmm_accession': str, 'E-value': float, 'score': float}
    df = pandas.read_csv(tsv, sep='\t', usecols=useful_column_names, dtype=column_types)
    df = df.loc[df['E-value'] <= max_evalue]
    df = df.sort_values(by='E-value', ascending=True, kind='mergesort').drop_duplicates(subset=['target_name']) 
    return df     


def hmmsearch(faa, hmm_library, stdout, tblout, domtblout, cpu, max_evalue):
    hmmsearch_args = ['hmmsearch', '--cpu', cpu, '--tblout', tblout, 
                      '--domtblout', domtblout,'-o', stdout, hmm_library, faa]
    subprocess.run(hmmsearch_args, check = True)
    
    header = ('target_name\ttgt_accession\thmm_name\thmm_accession\t'
              'E-value\tscore\tbias\t'
              'E-value_best_domain\tscore_best_domain\tbias_best_domain'
              '\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tdescription_of_target\n')
    
    nb_tabs = len(header.split('\t'))-1
    tblout_to_tsv(tblout, tblout[:-4] + '.tsv', header, nb_tabs)
    
    header = ('target_name\ttgt_accession\ttgt_len\thmm_name\thmm_accession\thmm_len\t'
              'fullseq_E-value\tfull_seq_score\tfull_seq_bias\t'
              'domain_nb\ttotal_nb_domains\tc-Evalue\ti-Evalue\tdomain_score\tdomain_bias\t'
              'hmm_start\thmm_end\ttgt_start\ttgt_end\tenv_start\tenv_end\tacc\tdescription_of_target\n')
    nb_tabs = len(header.split('\t'))-1
    tblout_to_tsv(domtblout, domtblout[:-4] + '.tsv', header, nb_tabs)
    
    results_df = filter_hmmsearch_results(tblout[:-4] + '.tsv', max_evalue)

    # return results_df[['target_name', 'hmm_name', 'hmm_accession', 'E-value', 'score']]
    return results_df[['target_name', 'hmm_name', 'hmm_accession', 'E-value', 'score']]


#######################################
# tprpred
#######################################
def filter_tprpred_results(tsv, min_proba):
    useful_column_names = ['Sequence', 'Score', 'P-value', 'Probab']
    column_types = {'Sequence': str, 'Score': float, 'P-value': float, 'Probab': str}
    df = pandas.read_csv(tsv, sep='\t', usecols=useful_column_names, dtype=column_types)
    df['Probab'] = df['Probab'].str[:-1].astype('float')
    df = df.loc[df['Probab'] >= min_proba]
    return df


def tprpred_to_hmmsearch_df_format(df):
    df = df[['Sequence', 'P-value', 'Score']]
    df.columns = ['target_name', 'E-value', 'score']
    # df['hmm_name'] = 'tpr2.8'
    # df['hmm_accession'] = 'Tprpred'
    # return df[['target_name', 'hmm_name', 'hmm_accession', 'E-value', 'score']]
    df['hmm_name'] = 'tprpred'
    return df[['target_name', 'hmm_name', 'E-value', 'score']]
    
    
def tprpred(faa, outfile, rrnpp_detector_dir, current_dir, min_proba):
    tprpred_dir = os.path.join(rrnpp_detector_dir, 'tprpred')
    os.chdir(tprpred_dir)
    tprpred_args = ['perl', 'tprpred.pl', faa, '-r', 'tpr2.8.pp']
    with open(outfile, mode = 'w') as f:
        subprocess.run(tprpred_args, check = True, stdout = f)
    f.close()
    os.chdir(current_dir)
    
    df = filter_tprpred_results(outfile, min_proba)
    df = tprpred_to_hmmsearch_df_format(df)
    return df


#######################################
# orfipy
#######################################
# def orfipy(fna, flanking_dict, parameters, cpu, working_dir):
#     min_len = parameters['propeptide_len_boundaries'][0] * 3
#     max_len = parameters['propeptide_len_boundaries'][1] * 3
#     orfipy_args = ['orfipy', '--outdir', working_dir, '--procs', cpu, '--ignore-case', '--dna', 'orfipy_micropeptides.fna', 
#                    '--pep', 'orfipy_micropeptides.faa', '--bed', 'orfipy_micropeptides.bed', '--strand', 'f', 
#                    '--min', str(min_len), '--max', str(max_len), '--start', ','.join(parameters['orfipy_start_codons']), fna]
#     subprocess.run(orfipy_args, check = True)
    
 
#######################################
# BlastP
#######################################
# apply thresholds and keep only the best homolog for each protein   
def filter_blastp_results(blastp_results, blastp_min_pident, blastp_min_pcover, blastp_max_evalue):
    if os.stat(blastp_results).st_size == 0:
        filtered_df = pandas.DataFrame(columns=['target_name', 'query_name', 'blastp_evalue', 'perc_identity', 'perc_cover'])
    else:
        df = pandas.read_csv(blastp_results, sep='\t', header=None)
        df.columns = ['query_name', 'target_name', 'qlen', 'slen', 'perc_identity', 'blastp_evalue', 'bitscore', 'qstart', 'qend', 'sstart', 'send']
        query_coverage = round((df['qend'] - df['qstart'] + 1) / df['qlen'] * 100, 1)
        target_coverage = round((df['send'] - df['sstart'] + 1) / df['slen'] * 100, 1)
        df['perc_cover'] = pandas.concat([query_coverage, target_coverage], axis=1).apply(min, axis=1)
        filtered_df = df.loc[(df['perc_identity'] >= blastp_min_pident) & (df['perc_cover'] >= blastp_min_pcover) & (df['blastp_evalue'] <= blastp_max_evalue)]
        filtered_df = filtered_df.sort_values(by = 'blastp_evalue', ascending = True, kind = 'mergesort').drop_duplicates(subset=['target_name'])
    return filtered_df[['target_name', 'query_name', 'blastp_evalue', 'perc_identity', 'perc_cover']]


def blastp(query_faa, target_faa, out_dir, cpu, blastp_min_pident, blastp_min_pcover, blastp_max_evalue):
    basename_target = os.path.basename(target_faa).split('.fa')[0]
    basename_query = os.path.basename(query_faa).split('.fa')[0]
    
    # makeblastdb
    database = os.path.join(out_dir, basename_target + '_blastpDB')
    makeblastdb_args = ['makeblastdb', '-dbtype', 'prot', '-in', target_faa, '-out', database]
    subprocess.run(makeblastdb_args, check = True)
    
    # blastp
    blastp_results = os.path.join(out_dir, basename_query + '_vs_' + basename_target + '_blastp_results.tsv')
    blastp_args = ['blastp', '-query', query_faa, '-db', database, '-out', blastp_results, '-num_threads', cpu, 
                   '-max_target_seqs', '1000000', '-outfmt', '6 qseqid sseqid qlen slen pident evalue bitscore qstart qend sstart send']
    subprocess.run(blastp_args, check = True)

    # clean blast
    df = filter_blastp_results(blastp_results, blastp_min_pident, blastp_min_pcover, blastp_max_evalue)
    return df
    

#######################################
# Signalp
#######################################
def signalp(faa, out_dir, rrnpp_detector_dir, current_dir):
    try:
        signalp_dir = os.path.dirname(shutil.which('signalp'))
    except:
        print('warning: the signalp binary is not found in your $PATH. RRNPP_detector will try to run built-in signalp version 5.0b Linux x86_64')
        signalp_dir = os.path.join(rrnpp_detector_dir, 'signalp-5.0b', 'bin')
    # signalp must be run from the directory where the binary is stored
    os.chdir(signalp_dir)
    signalp_args = ['./signalp', '-org', 'gram+', '-format', 'short', '-fasta', faa, '-prefix', os.path.join(out_dir, 'signalp')]
    subprocess.run(signalp_args, check = True)
    os.chdir(current_dir)
    df = pandas.read_csv(os.path.join(out_dir, 'signalp_summary.signalp5'), comment='#', sep='\t', header=None)
    df.columns = ['protein_id', 'Prediction', 'SP(Sec/SPI)', 'TAT(Tat/SPI)', 'LIPO(Sec/SPII)', 'OTHER', 'CS_Position']
    df['propeptide_score'] = df['SP(Sec/SPI)']
    return df[['protein_id', 'Prediction', 'propeptide_score', 'SP(Sec/SPI)', 'TAT(Tat/SPI)', 'LIPO(Sec/SPII)', 'OTHER', 'CS_Position']]

#######################################
# PrediSi
#######################################
def predisi(faa, out_dir, rrnpp_detector_dir, current_dir, min_likelihood):
    predisi_dir = os.path.join(rrnpp_detector_dir, 'predisi')
    os.chdir(predisi_dir)
    try:
        subprocess.call(['java', 'JSPP',
                         os.path.join(predisi_dir, 'matrices(SearchMatrix-objects)', 'gramp.smx'),
                         faa, os.path.join(out_dir, 'predisi_output.tsv')])
        df = pandas.read_csv(os.path.join(out_dir, 'predisi_output.tsv'), sep = "\t", header=None)
        df.columns = ['protein_id', 'CS_Position', 'Prediction', 'propeptide_score']
        df = df.loc[df['propeptide_score'] >= min_likelihood]
        # df = df.loc[df['Prediction'] == 'Y']
        df['Prediction'] = df['Prediction'].map({'Y': 'PrediSi_pos', 'N': 'PrediSi_neg_(but_high_score)'})
        df['protein_id'] = df['protein_id'].apply(lambda x: x.split(' ')[0])
        df['CS_Position'] = df['CS_Position'].apply(lambda x: "CS pos: " + str(x) + ':' + str(x + 1))
        return df
    except:
        print('Warning: Predisi exited with an error (only signalp results will be considered)')
        return None
    