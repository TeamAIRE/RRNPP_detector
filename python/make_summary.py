#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in genomes of bacteria and bacteriophages
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import pandas


def setdiff(long_list, short_list):
    s = set(short_list)
    only_in_long_list = [x for x in long_list if x not in s]
    return only_in_long_list
            

def add_coordinates_to_df(df, protein_dict, prefix, add_genomic_accession):
    proteins = df.iloc[:,0].tolist()
    start = list()
    end = list()
    strand = list()
    genomic_accession = list()
    for protein in proteins:
        start.append(protein_dict[protein]['start'])
        end.append(protein_dict[protein]['end'])
        strand.append(protein_dict[protein]['strand'])
        genomic_accession.append(protein_dict[protein]['genomic_accession'])
    df.insert(1, prefix + '_start', start)
    df.insert(2, prefix + '_end', end)
    df.insert(3, prefix + '_strand', strand)
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
    
    
def make_summary(hmm_tpr_df, hmm_tf_df, blast_df, signalp_df, protein_dict, position_dict):
    hmm_tpr_df = add_coordinates_to_df(hmm_tpr_df, protein_dict, 'r', True)
    signalp_df = add_coordinates_to_df(signalp_df, protein_dict, 'p', False)
    
    # receptor part
    hmm_tpr_df = hmm_tpr_df.rename(columns={'hmm_name': 'TPR_hmm_name', 'hmm_accession': 'TPR_hmm_accession', 'E-value': 'TPR_hmm_evalue'})
    hmm_tf_df = hmm_tf_df.rename(columns={'hmm_name': 'DNA-binding_hmm_name', 'hmm_accession': 'DNA-binding_hmm_accession', 'E-value': 'DNA-binding_hmm_evalue'})
    temp_df = pandas.merge(hmm_tpr_df, hmm_tf_df, how='left', left_on='target_name', right_on='target_name')
    receptors_df = pandas.merge(temp_df, blast_df, how='left', left_on='target_name', right_on='target_name')
    receptors_df = receptors_df.rename(columns={'target_name': 'receptor_id'})
    receptors_ids = receptors_df['receptor_id'].tolist()
    
    # intersect receptor part and propeptide part
    signalp_df = signalp_df.rename(columns={'protein_id': 'propeptide_id'})
    propeptides_ids = signalp_df['propeptide_id'].tolist()
    inter_df = pandas.DataFrame(columns = receptors_df.columns.tolist() + signalp_df.columns.tolist())
    for receptor in receptors_ids:
        genomic_accession = protein_dict[receptor]['genomic_accession']
        position = protein_dict[receptor]['position']
        flanking_orfs = list()
        if (genomic_accession, position-1) in position_dict:
            flanking_orfs.append(position_dict[(genomic_accession, position-1)]['protein_id'])
        if (genomic_accession, position+1) in position_dict:
            flanking_orfs.append(position_dict[(genomic_accession, position+1)]['protein_id'])
        for flanking_orf in flanking_orfs:
            if flanking_orf in propeptides_ids:
                suffix = signalp_df[signalp_df['propeptide_id'] == flanking_orf].reset_index(drop=True)
                prefix = receptors_df.loc[receptors_df['receptor_id'] == receptor].reset_index(drop=True)
                new_row = pandas.concat([prefix, suffix], axis=1)
                inter_df = inter_df.append(new_row)
    inter_df = inter_df.reset_index().drop('index', axis=1)
    # add intergenic distance and context columns
    final_df = add_intergenic_distance_and_context(inter_df)
    final_df = final_df.sort_values('genomic_accession', ascending = True, kind = 'mergesort')
    final_df = final_df.assign(qss_id=range(0, len(final_df)))
#    final_df.insert(0, 'genomic_accession', seq_id)
    return final_df
    
    