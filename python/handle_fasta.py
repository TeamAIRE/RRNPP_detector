#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in gram-positive bacteria and bacteriophages genomes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import os
import shutil
import sys


def filter_by_length(in_faa, out_faa, min_len, max_len):
    retained_proteins = list()
    with open(out_faa, mode='w') as outfile:
        l, header, sequence, protein_id = [0, '', '', '']
        with open(in_faa, mode='r') as infile:
            for line in infile:
                if line[0] == ">":
                    if l >= min_len and l <= max_len:
                        outfile.write(header)
                        outfile.write(sequence)
                        retained_proteins.append(protein_id)
                    header = line
                    protein_id = header[1:].split(' ', 1)[0]
                    sequence = ''
                    l = 0;
                else:
                    sequence = sequence + line
                    l = l + len(line.strip('\n'))
        infile.close()
        if l >= min_len and l <= max_len:
            outfile.write(header)
            outfile.write(sequence)
            retained_proteins.append(protein_id)
    outfile.close()
    return retained_proteins


def explicitize_duplicates(in_faa, out_dir, duplicates):
    list_protein_ids = list()
    header, sequence, protein_id = ['', '', '']
    header_splits = ['', '']
    out_faa = os.path.join(out_dir, 'tmp.faa')
    with open(out_faa, mode='w') as outfile:
        with open(in_faa, mode='r') as infile:
            for line in infile:
                if line[0] == ">":
                    # if a protein is visited and is in the list of query ids: write the sequence in the output faa
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


def filter_by_adjacency(propeptides_ids, receptors_ids, protein_dict, position_dict):
    if len(propeptides_ids) < len(receptors_ids):
        list_a = propeptides_ids
        list_b = receptors_ids
    else:
        list_a = receptors_ids
        list_b = propeptides_ids
    filt_list_a = list()
    filt_list_b = list()
    unique_orfs = dict()
    # get each protein of the shortest list
    for orf in list_a:
        genomic_accession = protein_dict[orf]['genomic_accession']
        position = protein_dict[orf]['position']
        # get its two flanking proteins
        flanking_orfs = list()
        if (genomic_accession, position-1) in position_dict:
            flanking_orfs.append(position_dict[(genomic_accession, position-1)]['protein_id'])
        if (genomic_accession, position+1) in position_dict:
            flanking_orfs.append(position_dict[(genomic_accession, position+1)]['protein_id'])
        # test if any flanking protein is in the longest list
        for flanking_orf in flanking_orfs:
            if flanking_orf in list_b:
                if not flanking_orf in unique_orfs:
                    filt_list_b.append(flanking_orf)
                    unique_orfs[flanking_orf] = True
                if not orf in unique_orfs:
                    filt_list_a.append(orf)
                    unique_orfs[orf] = True
    if len(propeptides_ids) < len(receptors_ids):
        return filt_list_a, filt_list_b
    else:
        return filt_list_b, filt_list_a
    

def subset_by_id(in_faa, out_faa, ids):
    k = 0
    # create a dictionary of visited proteins in order to skip reading fasta if all query proteins have been extracted from the faa 
    is_visited = dict()
    for protein_id in ids:
        is_visited[protein_id] = False
    with open(out_faa, mode='w') as outfile:
        header, sequence, protein_id = ['', '', '']
        with open(in_faa, mode='r') as infile:
            for line in infile:
                if line[0] == ">":
                    # if a protein is visited and is in the list of query ids: write the sequence in the output faa
                    if protein_id in ids:
                        outfile.write(header)
                        outfile.write(sequence)
                        k += 1
                        is_visited[protein_id] = True
                        if all(is_visited.values()):
                            break
                    header = line
                    protein_id = header[1:].split(' ', 1)[0]
                    sequence = ''
                else:
                    sequence = sequence + line
        infile.close()
        # handle last sequence of the input faa
        if protein_id in ids and not is_visited[protein_id]:
            outfile.write(header)
            outfile.write(sequence)
            k += 1
    outfile.close()
    if k == 0:
        sys.exit('error: unable to interesect annotation file and faa')
    elif k < len(ids):
        print('warning: could not intersect annotation file and faa for %d proteins' % (len(ids)-k))

