#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QSS_detector: a tool to detect quorum sensing systems
in bacterial genomes/metagenomes/contigs and mobile genetic elements
@author: Bernard Charles
"""

import pandas
            

def add_coordinates_to_df(df, protein_dict, prefix, add_genomic_accession):
    proteins = df.iloc[:,0].tolist()
    start = list()
    end = list()
    strand = list()
    genomic_accession = list()
    if prefix == 'p':
        rbs_bin = list()
        rbs_motif = list()
        rbs_spacer = list()
    for protein in proteins:
        start.append(protein_dict[protein]['start'])
        end.append(protein_dict[protein]['end'])
        strand.append(protein_dict[protein]['strand'])
        genomic_accession.append(protein_dict[protein]['genomic_accession'])
        if prefix == 'p':
            if 'RBS_bin' in protein_dict[protein]:
                rbs_bin.append(protein_dict[protein]['RBS_bin'])
                rbs_motif.append(protein_dict[protein]['RBS_motif'])
                rbs_spacer.append(protein_dict[protein]['RBS_spacer'])
            else:
                rbs_bin.append(None)
                rbs_motif.append(None)
                rbs_spacer.append(None)
    df.insert(1, prefix + '_start', start)
    df.insert(2, prefix + '_end', end)
    df.insert(3, prefix + '_strand', strand)
    if prefix == 'p':
        df.insert(4, prefix + '_RBS_bin', rbs_bin)
        df.insert(5, prefix + '_RBS_motif', rbs_motif)
        df.insert(6, prefix + '_RBS_spacer', rbs_spacer)
    if add_genomic_accession:
        df.insert(0, 'genomic_accession', genomic_accession)
    df = df.sort_values(by = prefix + '_start', ascending = True, kind = 'mergesort')
    return df


def add_intergenic_distance_and_context(df):
    intergenic_distances = list()
    contexts = list()
    for i in range(0, len(df)):
        r_start, r_end, r_strand, p_start, p_end, p_strand = df.loc[i, ['r_start', 'r_end', 'r_strand', 'p_start', 'p_end', 'p_strand']].tolist() 
        if r_start < p_start:
            intergenic_distances.append(p_start - r_end - 1)
        else:
            intergenic_distances.append(r_start - p_end - 1)
        if r_strand == p_strand:
            if r_strand == '+':
                if r_start < p_start:
                    contexts.append('(====> =>)')
                else:
                    contexts.append('(=> ====>)')
            elif r_strand == '-':
                if r_start < p_start:
                    contexts.append('(<==== <=)')
                else:
                    contexts.append('(<= <====)')
            else:
                contexts.append('')
        else:
            if r_strand == '+' and p_strand == '-':
                if r_start < p_start:
                    contexts.append('(====> <=)')
                else:
                    contexts.append('(<= ====>)')
            elif r_strand == '-' and p_strand == '+':
                if r_start < p_start:
                    contexts.append('(<==== =>)')
                else:
                    contexts.append('(=> <====)')
            else:
                contexts.append('')
    df.insert(1, 'genomic_context', contexts)
    df.insert(2, 'intergenic_distance', intergenic_distances)
    return df


def merge_signalp_predisi(signalp_df, predisi_df):
    idx = signalp_df[signalp_df['protein_id'].isin(predisi_df['protein_id'])].index
    signalp_df.loc[idx,['CS_Position', 'Prediction', 'propeptide_score']] = predisi_df[['CS_Position', 'Prediction', 'propeptide_score']].values
    return signalp_df
    
    
    
def make_summary(tpr_df, tf_df, blast_df, signalp_df, protein_dict, cognate_dict):
    tpr_df = add_coordinates_to_df(tpr_df, protein_dict, 'r', True)
    signalp_df = add_coordinates_to_df(signalp_df, protein_dict, 'p', False)
    signalp_df.to_csv('propeptides_summary.tsv', sep='\t', header=True, index=False, mode='w')
    
    # receptor part
    tpr_df = tpr_df.rename(columns={'hmm_name': 'TPR_hmm_name', 'hmm_accession': 'TPR_hmm_accession', 'E-value': 'TPR_hmm_evalue'})
    tf_df = tf_df.rename(columns={'hmm_name': 'DNA-bd_hmm_name', 'hmm_accession': 'DNA-bd_hmm_accession', 'E-value': 'DNA-bd_hmm_evalue'})
    temp_df = pandas.merge(tpr_df.drop('score', axis=1), tf_df.drop('score', axis=1), how='left', left_on='target_name', right_on='target_name')
    receptors_df = pandas.merge(temp_df, blast_df, how='left', left_on='target_name', right_on='target_name')
    receptors_df = receptors_df.rename(columns={'target_name': 'receptor_id'})
    receptors_ids = receptors_df['receptor_id'].tolist()
    
    # intersect receptor part and propeptide part
    signalp_df = signalp_df.rename(columns={'protein_id': 'propeptide_id'})
    signalp_df.to_csv('signalp_summary.tsv', sep='\t', header=True, index=False, mode='w')
    inter_df = pandas.DataFrame(columns = receptors_df.columns.tolist() + signalp_df.columns.tolist())
    
    for receptor in receptors_ids:
        if receptor in cognate_dict:
            for propeptide in cognate_dict[receptor]['cognate']:
                suffix = signalp_df[signalp_df['propeptide_id'] == propeptide].reset_index(drop=True)
                prefix = receptors_df.loc[receptors_df['receptor_id'] == receptor].reset_index(drop=True)
                new_row = pandas.concat([prefix, suffix], axis=1)
                inter_df = pandas.concat([inter_df, new_row], axis=0)
    inter_df = inter_df.reset_index().drop('index', axis=1)
    # add intergenic distance and context columns
    final_df = add_intergenic_distance_and_context(inter_df)
    final_df = final_df.sort_values('genomic_accession', ascending = True, kind = 'mergesort')
#    final_df.insert(0, 'genomic_accession', seq_id)
    return final_df


def make_iterative_search_summary(receptors_blast_df, propeptides_blast_df, protein_dict, cognate_dict):
    receptors_blast_df = add_coordinates_to_df(receptors_blast_df, protein_dict, 'r', True)
    propeptides_blast_df = add_coordinates_to_df(propeptides_blast_df, protein_dict, 'p', True)
    propeptides_blast_df = propeptides_blast_df.drop('genomic_accession', axis = 1)
    receptors_blast_df = receptors_blast_df.rename(columns={'target_name': 'homologous_receptor', 'query_name': 'query_receptor', 
                                                            'blastp_evalue': 'r_blastp_evalue', 'perc_identity': 'r_perc_identity', 
                                                            'perc_cover': 'r_perc_cover'})
    propeptides_blast_df = propeptides_blast_df.rename(columns={'target_name': 'homologous_propeptide', 'query_name': 'query_propeptide', 
                                                            'blastp_evalue': 'p_blastp_evalue', 'perc_identity': 'p_perc_identity', 
                                                            'perc_cover': 'p_perc_cover'})
    inter_df = pandas.DataFrame(columns = receptors_blast_df.columns.tolist() + propeptides_blast_df.columns.tolist())
    
    receptors_ids = receptors_blast_df['homologous_receptor'].tolist()
    propeptides_dict = dict.fromkeys(propeptides_blast_df['homologous_propeptide'].tolist(), True)
    
    for receptor in receptors_ids:
        if receptor in cognate_dict:
            for propeptide in cognate_dict[receptor]['cognate']:
                if propeptide in propeptides_dict:
                    suffix = propeptides_blast_df[propeptides_blast_df['homologous_propeptide'] == propeptide].reset_index(drop=True)
                    prefix = receptors_blast_df.loc[receptors_blast_df['homologous_receptor'] == receptor].reset_index(drop=True)
                    new_row = pandas.concat([prefix, suffix], axis=1)
                    inter_df = pandas.concat([inter_df, new_row], axis=0)
    inter_df = inter_df.reset_index().drop('index', axis=1)
    # add intergenic distance and context columns
    final_df = add_intergenic_distance_and_context(inter_df)
    final_df = final_df.sort_values('genomic_accession', ascending = True, kind = 'mergesort')
    return final_df