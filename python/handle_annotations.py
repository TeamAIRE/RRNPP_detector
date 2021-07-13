#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in genomes of bacteria and bacteriophages
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import itertools
import os
import pandas
import subprocess
import sys


def prodigal(fna, out_dir):
    # produce a faa and a gff from a fna
    gff = os.path.join(out_dir, 'predicted_ORFs.gff')
    faa = os.path.join(out_dir, 'predicted_proteins.faa')
    prodigal_args = ['prodigal', '-i', fna, '-f', 'gff', '-o', gff, '-a', faa]
    subprocess.run(prodigal_args, check = True)
    print('  Prodigal done!')
    return gff, faa


def check_anno_coordinates(df):
    if not df['start'].dtype == 'int':
        sys.exit('error: column \'start\' (field n°%d) in the input annotation file is not of type integer' % df.columns.get_loc('start'))
    if not df['end'].dtype == 'int':
        sys.exit('error: column \'end\' (field n°%d) in the input annotation file is not of type integer' % df.columns.get_loc('end'))
    if not df['strand'].str.match(r'^[+-.]$').all():
       sys.exit('error: column \'strand\' (field n°%d) in the input annotation file is not filled only with the expected \'+\', \'-\' or \'.\' characters' % df.columns.get_loc('strand'))
    

def feature_tbl_to_dict(feature_tbl):
    ft_column_names = ['feature', 'class', 'assembly', 'assembly_unit', 'seq_type', 'chromosome', 'genomic_accession', 'start', 'end', 'strand', 'protein_id',
                       'non-redundant_refseq', 'related_accession','name', 'symbol', 'GeneID', 'locus_tag', 'feature_interval_length', 'product_length', 'attributes']
    useful_column_names = ['feature', 'class', 'genomic_accession', 'start', 'end', 'strand', 'protein_id']
    column_types = {'feature': str, 'class': str, 'genomic_accession': str, 'start': int, 'end': int, 'strand': str, 'protein_id': str}
    try:
        df = pandas.read_csv(feature_tbl, sep='\t', comment='#', header=None, names=ft_column_names, usecols=useful_column_names, dtype=column_types)
    except Exception as e:
        print('EXIT: the annotation file does not comply with the feature table format requirements.')
        print('Reading the feature table exited with the following error:')
        sys.exit(e)
    if not df['strand'].str.match(r'^[+-.]$').all():
       sys.exit('error: column \'strand\' (field n°%d) in the input annotation file is not filled only with the expected \'+\', \'-\' or \'.\' characters' % df.columns.get_loc('strand'))
    # retain only ORFs that encode proteins
    orf_df = df.loc[(df['feature'] == 'CDS') & (df['class'] == 'with_protein')].reset_index()
    # check there are cds in the feature table
    if orf_df.empty:
        sys.exit('error: columns \'features\' and \'class\' (fields n°1 and 2) of the input feature table never correspond to the \'CDS\twith_protein\' string')
    # drop rows where protein_id is null (usually correspond to pseudogenes)
    orf_df = orf_df.dropna(subset=['protein_id'])
    # make sure table is sorted
    orf_df = orf_df.sort_values(['genomic_accession', 'start'], ascending = True, kind = 'mergesort')
    # make sure all protein ids are unique (add suffix __dupli_nb1, __dupli_nb2, __dupli_nb3 ... otherwise)
    dupli_dict = dict()
    if not orf_df['protein_id'].is_unique:
        occurences = orf_df['protein_id'].value_counts()
        duplicates = occurences[occurences > 1]
        duplicated_proteins = duplicates.keys()
        for protein in duplicated_proteins:
            new_values = [protein + '__dupli_nb' + str(i+1) for i in range(duplicates[protein])]
            dupli_dict[protein] = new_values
            orf_df.loc[(orf_df['protein_id'] == protein), 'protein_id'] = new_values
    # create a position column based on the order of the ORFs
    orf_df = orf_df.assign(position=range(0, len(orf_df)))
    # create dictionaries from the data frame
    position_dict = orf_df[['genomic_accession', 'position', 'protein_id']].set_index(['genomic_accession', 'position']).T.to_dict()
    protein_dict = orf_df[['genomic_accession', 'protein_id', 'position', 'start', 'end', 'strand']].set_index('protein_id').T.to_dict()
    return dupli_dict, position_dict, protein_dict
    

def gff_attributes_to_columns(df):
    attribute_df = pandas.DataFrame(df['attributes'].apply(
        lambda attributes: dict([key_value_pair.split(sep='=', maxsplit=1) for key_value_pair in attributes.strip(';').split(';')])))
    attribute_df.columns = ['at_dic']
    attribute_df['at_dic_keys'] = attribute_df['at_dic'].apply(lambda at_dic: list(at_dic.keys()))
    merged_attribute_list = list(itertools.chain.from_iterable(attribute_df['at_dic_keys']))
    nonredundant_list = sorted(list(set(merged_attribute_list)))
    final_df = df.loc[:, "genomic_accession":"strand"]
    for atr in nonredundant_list:
        final_df[atr] = attribute_df["at_dic"].apply(lambda at_dic: at_dic.get(atr))
    return final_df
    

def gff_to_dict(gff):
    gff_column_names = ['genomic_accession', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    useful_column_names = ['genomic_accession', 'type', 'start', 'end', 'strand', 'attributes']
    column_types = {'genomic_accession': str, 'type': str, 'start': int, 'end': int, 'strand': str, 'attributes': str}
    try:
        df = pandas.read_csv(gff, sep='\t', comment='#', header=None, names=gff_column_names, usecols=useful_column_names, dtype=column_types)
    except Exception as e:
        print('EXIT: the annotation file does not comply with gff format requirements.')
        print('Reading the gff exited with the following error:')
        sys.exit(e)
    if not df['strand'].str.match(r'^[+-.]$').all():
       sys.exit('error: column \'strand\' (field n°%d) in the input annotation file is not filled only with the expected \'+\', \'-\' or \'.\' characters' % df.columns.get_loc('strand'))
    # retain only ORFs that encode proteins
    orf_df = df.loc[df['type'] == 'CDS']
    # check there are cds in the gff
    if orf_df.empty:
        sys.exit('error: column \'type\' (fields n°3) of the input gff never corresponds to the \'CDS\' string')
    # turn attributes string into columns
    orf_df = gff_attributes_to_columns(orf_df)
    # check a protein_id or a id column is present in the data frame
    if not 'protein_id' in orf_df:
        if 'ID' in orf_df:
            temp_id = orf_df['ID']
        elif 'id' in orf_df:
            temp_id = orf_df['id']
        else:
            sys.exit('Neither the tag \'protein_id\' nor \'ID\' are present in the column \'attributes\' of the input gff')
        orf_df['protein_id'] = orf_df['genomic_accession'] + temp_id.apply(lambda x: '_' + x.split('_')[1])  
    # drop rows where protein_id is null (usually correspond to pseudogenes)
    orf_df = orf_df.dropna(subset=['protein_id'])
    # make sure table is sorted
    orf_df = orf_df.sort_values(['genomic_accession', 'start'], ascending = True, kind = 'mergesort')
    # make sure all protein ids are unique (add suffix __dupli_nb1, __dupli_nb2, __dupli_nb3 ... otherwise)
    dupli_dict = dict()
    if not orf_df['protein_id'].is_unique:
        occurences = orf_df['protein_id'].value_counts()
        duplicates = occurences[occurences > 1]
        duplicated_proteins = duplicates.keys()
        for protein in duplicated_proteins:
            new_values = [protein + '__dupli_nb' + str(i+1) for i in range(duplicates[protein])]
            orf_df.loc[(orf_df['protein_id'] == protein), 'protein_id'] = new_values
            dupli_dict[protein] = new_values
    # create a position column based on the order of the ORFs
    orf_df = orf_df.assign(position=range(0, len(orf_df)))
    # create dictionaries from the data frame
    position_dict = orf_df[['genomic_accession', 'position', 'protein_id']].set_index(['genomic_accession', 'position']).T.to_dict()
    protein_dict = orf_df[['genomic_accession', 'protein_id', 'position', 'start', 'end', 'strand']].set_index('protein_id').T.to_dict()
    return dupli_dict, position_dict, protein_dict


def make_small_annotation_file(annofile, test_file):
    max_lines = 50
    l = 0
    with open(test_file, mode = 'w') as outfile:
         with open(annofile, mode='r') as infile:
             while l < max_lines:
                 line = infile.readline()
                 if not line[0] == '#':
                     l += 1
                     outfile.write(line)
         infile.close()
    outfile.close()


# the two following test functions will be called only with a chunk annotation file and will return precise errors regarding gff or ft format 
def test_feature_tbl(feature_tbl):
    # read table
    df = pandas.read_csv(feature_tbl, sep='\t', comment='#', header=None)                  
    # check number of columns is 20                    
    if len(df.columns) != 20:
        sys.exit('EXIT: Feature table \'%s\' is expected to have 20 columns. %d found' % (feature_tbl, len(df.columns)))   
    # make header
    df.columns = ['feature', 'class', 'assembly', 'assembly_unit', 'seq_type', 'chromosome', 
                  'genomic_accession', 'start', 'end', 'strand', 'protein_id', 
                  'non-redundant_refseq', 'related_accession','name', 'symbol', 'GeneID', 
                  'locus_tag', 'feature_interval_length', 'product_length', 'attributes']
    # some quality check
    check_anno_coordinates(df)
    # retain only ORFs that encode proteins
    orf_df = df.loc[(df['feature'] == 'CDS') & (df['class'] == 'with_protein')].reset_index()
    # check there are cds in the feature table
    if orf_df.empty:
        sys.exit('error: columns \'features\' and \'class\' (fields n°1 and 2) of the input feature table never correspond to the \'CDS\twith_protein\' string')

    
def test_gff(gff):
    # read table
    df = pandas.read_csv(gff, sep='\t', comment='#', header=None)  
    # check number of columns is 20                    
    if len(df.columns) != 9:
        sys.exit('EXIT: Gff \'%s\' is expected to have 9 columns. %d found' % (gff, len(df.columns)))
    # make header
    df.columns = ['genomic_accession', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    # some quality check
    check_anno_coordinates(df)
    # retain only ORFs that encode proteins
    orf_df = df.loc[df['type'] == 'CDS']
    # check there are cds in the gff
    if orf_df.empty:
        sys.exit('error: column \'type\' (fields n°3) of the input gff never corresponds to the \'CDS\' string')
    # turn attributes string into columns
    orf_df = gff_attributes_to_columns(orf_df)
    # check a protein_id or a id column is present in the data frame
    if not 'protein_id' in orf_df:
        if 'ID' in orf_df:
            temp_id = orf_df['ID']
        elif 'id' in orf_df:
            temp_id = orf_df['id']
        else:
            sys.exit('Neither the tag \'protein_id\' nor \'ID\' are present in the column \'attributes\' of the input gff')
        orf_df['protein_id'] = orf_df['genomic_accession'] + temp_id.apply(lambda x: '_' + x.split('_')[1])
    

def test_annotation_file(annotation_file, annotation_format, out_dir):
    # the gff_to_dict or feature_tbl_to_dict function will output error messages if annofile doesn't comply with gff or ft requirements
    if annotation_format == 'gff':
        test_file = os.path.join(out_dir, "small_test_gff.gff")
        make_small_annotation_file(annotation_file, test_file)
        test_gff(test_file)
    elif annotation_format == 'feature_tbl':
        test_file = os.path.join(out_dir, "small_test_feature_table.txt")
        make_small_annotation_file(annotation_file, test_file)
        test_feature_tbl(test_file)
        

def split_annotation_file(annotation_file, annotation_format, out_dir):
    if annotation_format == 'gff':
        extension = '.gff'
    elif annotation_format == 'feature_tbl':
        extension = '.txt'
    old_seq_id = ''
    list_seq_ids = list()
    with open(annotation_file, mode = 'r') as infile:
        for line in infile:
            if line[0] != '#':
                if annotation_format == 'gff':
                    seq_id = line.split('\t', 1)[0]
                elif annotation_format == 'feature_tbl':
                    seq_id = line.split('\t', 7)[6]
                if seq_id != old_seq_id:
                    list_seq_ids.append(seq_id)
                    if old_seq_id != '':
                        outfile.close()
                    outfile = open(os.path.join(out_dir, seq_id + extension), mode = 'w')
                    old_seq_id = seq_id
                outfile.write(line)
    infile.close()
    outfile.close()
    return list_seq_ids
   
                         