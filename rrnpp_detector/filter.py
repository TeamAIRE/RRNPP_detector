#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-type quorum sensing systems in chromosomes, plasmids and phages of Firmicutes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import sys


def filter_by_length(in_faa, out_faa, boundaries, stop_codon = False):
    if stop_codon:
        min_len = boundaries[0] + 1
        max_len = boundaries[1] + 1
    else:
        min_len = boundaries[0] + 1
        max_len = boundaries[1] + 1
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


def filter_by_adjacency(proteins_a, proteins_b, nb_intercalar_genes, protein_dict, position_dict):
    if len(proteins_a) <= len(proteins_b):
        ref_list = proteins_a
        vs_dict = dict.fromkeys(proteins_b, True)
    else:
        ref_list = proteins_b
        vs_dict = dict.fromkeys(proteins_a, True)
    filt_ref_list = list()
    filt_vs_list = list()
    couple_dict = dict() # accessed by the id 'c', return a list with the two components forming a couple
    cognate_dict = dict() #accessed by a protein id, return the list of cognate components
    # (because can be flanked on both side by a component) as well as the couple key
    c = 0
    # get each protein of the shortest list
    for orf in ref_list:
        if orf in protein_dict:
            genomic_accession = protein_dict[orf]['genomic_accession']
            position = protein_dict[orf]['position']
            # get its flanking proteins
            flanking_orfs = list()
            for i in range(1, nb_intercalar_genes+2):
                if (genomic_accession, position-i) in position_dict:
                    flanking_orfs.append(position_dict[(genomic_accession, position-i)]['protein_id'])
                if (genomic_accession, position+i) in position_dict:
                    flanking_orfs.append(position_dict[(genomic_accession, position+i)]['protein_id'])
            # test if any flanking protein is in the longest list
            for flanking_orf in flanking_orfs:
                if flanking_orf in vs_dict:
                    c += 1
                    couple_dict[c] = [orf, flanking_orf]
                    if not flanking_orf in cognate_dict:
                        cognate_dict[flanking_orf] = dict()
                        cognate_dict[flanking_orf]['cognate'] = list()
                        cognate_dict[flanking_orf]['couple_id'] = list()
                        filt_vs_list.append(flanking_orf)
                    cognate_dict[flanking_orf]['cognate'].append(orf)
                    cognate_dict[flanking_orf]['couple_id'].append(c)
                    if not orf in cognate_dict:
                        cognate_dict[orf] = dict()
                        cognate_dict[orf]['cognate'] = list()
                        cognate_dict[orf]['couple_id'] = list()
                        filt_ref_list.append(orf)
                    cognate_dict[orf]['cognate'].append(flanking_orf)
                    cognate_dict[orf]['couple_id'].append(c)
        else:
            print('  warning: %s not in protein_dict' % orf)
    if len(proteins_a) <= len(proteins_b):
        return filt_ref_list, filt_vs_list, couple_dict, cognate_dict
    else:
        return filt_vs_list, filt_ref_list, couple_dict, cognate_dict
    

def subset_by_id(in_faa, subset_faa, protein_ids, junk_faa = None):
    # if filtered_faa is provided, then the filtered proteins will be written in this file
    # for instance, this is useful for LuxR solos 
    ids = dict.fromkeys(protein_ids, True)
    n = len(ids)
    k = 0
    # when k == n, all protein_ids to retain have been visited -> skip reading fasta
    with open(subset_faa, mode='w') as subset_file:
        if junk_faa:
            junk_file = open(junk_faa, mode='w')
        header, sequence, protein_id = ['', '', '']
        with open(in_faa, mode='r') as infile:
            for line in infile:
                if line[0] == ">":
                    # if a protein is visited and is in the list of query ids: write the sequence in the output faa
                    if protein_id in ids:
                        k += 1
                        if k == n and not junk_faa:
                            break
                        subset_file.write(header)
                        subset_file.write(sequence)
                    elif protein_id != '' and junk_faa:
                        junk_file.write(header)
                        junk_file.write(sequence)
                    header = line
                    protein_id = header[1:].split(' ', 1)[0]
                    sequence = '' 
                else:
                    sequence = sequence + line
        infile.close()
        # handle last sequence of the input faa
        if protein_id in ids:
            subset_file.write(header)
            subset_file.write(sequence)
            k += 1
        elif junk_faa:
            junk_file.write(header)
            junk_file.write(sequence)
    subset_file.close()
    if junk_faa:
        junk_file.close()
    if k == 0:
        sys.exit('error: unable to intersect annotation file and faa')
    elif k < n:
        print('warning: could not intersect annotation file and faa for %d proteins' % (n-k))
        

def get_only_partners(anchor_proteins, cognate_dict):
    reference_proteins = list()
    cognate_proteins = list()
    for p in anchor_proteins:
        if p in cognate_dict:
            reference_proteins.append(p)
            cognate_proteins += cognate_dict[p]['cognate']
    return reference_proteins, list(set(cognate_proteins))


def setdiff(long_list, short_list):
    s = set(short_list)
    only_in_long_list = [x for x in long_list if x not in s]
    return only_in_long_list

