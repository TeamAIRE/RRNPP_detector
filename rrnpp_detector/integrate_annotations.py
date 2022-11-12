#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-type quorum sensing systems in chromosomes, plasmids and phages of Firmicutes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import itertools
import os
import pandas
import shutil
import sys
from rrnpp_detector.wrappers import prodigal


def explicitize_duplicates(in_faa, out_dir, duplicates):
    list_protein_ids = list()
    header, sequence, protein_id = ['', '', '']
    header_splits = ['', '']
    out_faa = os.path.join(out_dir, 'tmp.faa')
    with open(out_faa, mode='w') as outfile:
        with open(in_faa, mode='r') as infile:
            for line in infile:
                if line[0] == ">":
                    # if a protein is visited and is in the list of query ids: 
                    # write the sequence in the output faa
                    if protein_id != '':
                        if protein_id in duplicates:
                            for unique_id in duplicates[protein_id]:
                                header = '>' + unique_id + ' ' + header_splits[1]
                                outfile.write(header)
                                outfile.write(sequence)
                                list_protein_ids.append(unique_id)
                        else:
                            outfile.write(header)
                            outfile.write(sequence)
                            list_protein_ids.append(protein_id)
                    header = line
                    header_splits = header[1:].split(' ', 1)
                    protein_id = header_splits[0]
                    sequence = ''
                else:
                    sequence = sequence + line
        infile.close()
        # handle last sequence of the input faa
        if protein_id in duplicates:
            for unique_id in duplicates[protein_id]:
                header = '>' + unique_id + ' ' + header_splits[1]
                outfile.write(header)
                outfile.write(sequence)
                list_protein_ids.append(unique_id)
        else:
            outfile.write(header)
            outfile.write(sequence)
            list_protein_ids.append(protein_id)
    outfile.close()
    shutil.move(out_faa, in_faa)
    return list_protein_ids


def check_anno_coordinates(df):
    if not df['start'].dtype == 'int':
        sys.exit('error: column \'start\' (field n°%d) '
                 'in the input annotation file is not of type integer' 
                 % df.columns.get_loc('start'))
    if not df['end'].dtype == 'int':
        sys.exit('error: column \'end\' (field n°%d) '
                 'in the input annotation file is not of type integer' 
                 % df.columns.get_loc('end'))
    if not df['strand'].str.match(r'^[+-.]$').all():
       sys.exit('error: column \'strand\' (field n°%d) '
                'in the input annotation file is not filled '
                'only with the expected \'+\', \'-\' or \'.\' characters'
                % df.columns.get_loc('strand'))
    

def feature_tbl_to_dict(feature_tbl):
    ft_column_names = ['feature', 'class', 'assembly', 'assembly_unit', 
                       'seq_type', 'chromosome', 'genomic_accession', 
                       'start', 'end', 'strand', 'protein_id',
                       'non-redundant_refseq', 'related_accession',
                       'name', 'symbol', 'GeneID', 'locus_tag', 
                       'feature_interval_length', 'product_length', 
                       'attributes']
    useful_column_names = ['feature', 'class', 'genomic_accession', 
                           'start', 'end', 'strand', 'protein_id']
    column_types = {'feature': str, 'class': str, 'genomic_accession': str, 
                    'start': int, 'end': int, 'strand': str, 'protein_id': str}
    try:
        df = pandas.read_csv(feature_tbl, sep='\t', comment='#', 
                             header=None, names=ft_column_names, 
                             usecols=useful_column_names, dtype=column_types)
    except Exception as e:
        print('EXIT: the annotation file does not comply '
              'with the feature table format requirements.')
        print('Reading the feature table exited with the following error:')
        sys.exit(e)
    if not df['strand'].str.match(r'^[+-.?]$').all():
       sys.exit('error: column \'strand\' (field n°%d) '
                'in the input annotation file is not filled '
                'only with the expected \'+\', \'-\' or \'.\' characters' 
                % df.columns.get_loc('strand'))
    # retain only ORFs that encode proteins
    orf_df = df.loc[(df['feature'] == 'CDS') & (df['class'] == 'with_protein')].reset_index()
    # check there are cds in the feature table
    if orf_df.empty:
        sys.exit('error: columns \'features\' and \'class\' '
                 '(fields n°1 and 2) of the input feature table '
                 'never correspond to the \'CDS\twith_protein\' string')
    # drop rows where protein_id is null (usually correspond to pseudogenes)
    orf_df = orf_df.dropna(subset=['protein_id'])
    # make sure table is sorted
    # orf_df = orf_df.sort_values(['genomic_accession', 'start'], ascending = True, kind = 'mergesort')
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
    

def gff_attributes_to_columns(df, max_attrs = None):
    attribute_df = pandas.DataFrame(df['attributes'].apply(
        lambda attributes: dict([key_value_pair.split(sep='=', maxsplit=1) for key_value_pair in attributes.strip(';').split(';')])))
    attribute_df.columns = ['at_dic']
    attribute_df['at_dic_keys'] = attribute_df['at_dic'].apply(lambda at_dic: list(at_dic.keys()))
    merged_attribute_list = list(itertools.chain.from_iterable(attribute_df['at_dic_keys']))
    nonredundant_list = sorted(list(set(merged_attribute_list)))
    df = df.loc[:, "genomic_accession":"strand"]
    for atr in nonredundant_list:
        df[atr] = attribute_df["at_dic"].apply(lambda at_dic: at_dic.get(atr))
    return df
    

def gff_to_dict(gff):
    commented_lines = list()
    k = 0
    with open(gff, mode='r') as f:
        for line in f:
            if line[0] == "#":
                commented_lines.append(k)
            k += 1
    f.close()
    
    gff_column_names = ['genomic_accession', 'source', 'type', 'start', 
                        'end', 'score', 'strand', 'phase', 'attributes']
    useful_column_names = ['genomic_accession', 'type', 'start', 
                           'end', 'strand', 'attributes']
    column_types = {'genomic_accession': str, 'type': str, 'start': int, 
                    'end': int, 'strand': str, 'attributes': str}
            
    try:
        df = pandas.read_csv(gff, sep='\t', skip_blank_lines=True, skiprows=commented_lines, header=None, names=gff_column_names, usecols=useful_column_names, dtype=column_types)
    except Exception as e:
        print('EXIT: the annotation file does not comply '
              'with gff format requirements.')
        print('Reading the gff exited with the following error:')
        sys.exit(e)
    if not df['strand'].str.match(r'^[+-.?]$').all():
       sys.exit('error: column \'strand\' (field n°%d) '
                'in the input annotation file is not filled '
                'only with the expected \'+\', \'-\' or \'.\' characters' 
                % df.columns.get_loc('strand'))
       
    # retain only ORFs that encode proteins
    df = df.loc[df['type'] == 'CDS']
    
    # check there are cds in the gff
    if df.empty:
        sys.exit('error: column \'type\' (fields n°3) '
                 'of the input gff never corresponds to the \'CDS\' string')
        
    # turn attributes string into columns
    # df = gff_attributes_to_columns(df)
    
    
    if (df['attributes'].iat[0] + ";protein_id=").split(sep = "protein_id=", maxsplit=1)[1].split(";")[0]:
        df['protein_id'] = df['attributes'].apply(
            lambda attributes: (attributes + ";protein_id=").split(sep = "protein_id=", maxsplit=1)[1].split(";")[0])
    elif (df['attributes'].iat[1] + ";Name=").split(sep = "Name=", maxsplit=1)[1].split(";")[0]:
        df['protein_id'] = df['attributes'].apply(
            lambda attributes: (attributes + ";Name=").split(sep = "Name=", maxsplit=1)[1].split(";")[0])
    elif (df['attributes'].iat[0] + ";ID=").split(sep = "ID=", maxsplit=1)[1].split(";")[0]:
        df['protein_id'] = df['genomic_accession'] + df['attributes'].apply(
            lambda attributes: "_" + (attributes + ";ID=").split(sep = "ID=", maxsplit=1)[1].split(";")[0].split('_')[1])
    elif (df['attributes'].iat[0] + ";id=").split(sep = "id=", maxsplit=1)[1].split(";")[0]:
        df['protein_id'] = df['genomic_accession'] +  df['attributes'].apply(
            lambda attributes: "_" + (attributes + ";id=").split(sep = "id=", maxsplit=1)[1].split(";")[0].split('_')[1])
    else:
        sys.exit('Neither the tag \'protein_id\' nor \'ID\' are present in the column \'attributes\' of the input gff')
            
    if not 'protein_id' in df:
        sys.exit('the tag \'protein_id\ is absent from the column \'attributes\' of the input gff for at least one CDS entry')
    
    # drop rows where protein_id is null (usually correspond to pseudogenes)
    df = df.dropna(subset=['protein_id'])
    
    # make sure table is sorted
    # orf_df = orf_df.sort_values(['genomic_accession', 'start'], ascending = True, kind = 'mergesort')
    # make sure all protein ids are unique (add suffix __dupli_nb1, __dupli_nb2, __dupli_nb3 ... otherwise)
    dupli_dict = dict()
    if not df['protein_id'].is_unique:
        occurences = df['protein_id'].value_counts()
        duplicates = occurences[occurences > 1]
        duplicated_proteins = duplicates.keys()
        for protein in duplicated_proteins:
            new_values = [protein + '__dupli_nb' + str(i+1) for i in range(duplicates[protein])]
            df.loc[(df['protein_id'] == protein), 'protein_id'] = new_values
            dupli_dict[protein] = new_values
            
    # create a position column based on the order of the ORFs
    df = df.assign(position=range(0, len(df)))
    df = df[useful_column_names[:-1] + ['protein_id', 'position']]
    
    # create dictionaries from the data frame
    position_dict = df[['genomic_accession', 'position', 'protein_id']].set_index(['genomic_accession', 'position']).T.to_dict()
    protein_dict = df[['genomic_accession', 'protein_id', 'position', 'start', 'end', 'strand']].set_index('protein_id').T.to_dict()
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
    df.columns = ['genomic_accession', 'source', 'type', 'start', 'end', 
                  'score', 'strand', 'phase', 'attributes']
    # some quality check
    check_anno_coordinates(df)
    # retain only ORFs that encode proteins
    orf_df = df.loc[df['type'] == 'CDS']
    # check there are cds in the gff
    if orf_df.empty:
        sys.exit('error: column \'type\' (fields n°3) '
                 'of the input gff never corresponds to the \'CDS\' string')
    # turn attributes string into columns
    orf_df = gff_attributes_to_columns(orf_df)
    # check a protein_id or a id column is present in the data frame
    if not 'protein_id' in orf_df:
        if 'ID' in orf_df:
            temp_id = orf_df['ID']
        elif 'id' in orf_df:
            temp_id = orf_df['id']
        else:
            sys.exit('Neither the tag \'protein_id\' nor \'ID\' '
                     'are present in the column \'attributes\' of the input gff')
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
        extension = '.tsv'
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


def get_annotation_config(faa, fna, out_dir, args):
    if faa:
        test_dir = os.path.join(out_dir, "tmp_annotation_directory")
        os.makedirs(test_dir)
        # Check annotation file complies with requirements of its format
        if args.feature_tbl:
            annotation_format = 'feature_tbl'
            annotation_file = args.feature_tbl
        elif args.gff:
            annotation_format = 'gff'
            annotation_file = args.gff
        test_annotation_file(annotation_file, annotation_format, test_dir)
        shutil.rmtree(test_dir)
    else:
        prodigal_dir = os.path.join(out_dir, "annotations")
        os.makedirs(prodigal_dir)
        print('-------------------------------------')
        print('* Running Prodigal to detect ORFs in the target genome(s) ...')
        annotation_file, faa = prodigal(fna, prodigal_dir)
        annotation_format = 'gff'
    return faa, annotation_file, annotation_format
   
                         
