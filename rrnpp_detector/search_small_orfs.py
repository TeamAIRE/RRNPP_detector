#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-type quorum sensing systems in chromosomes, plasmids and phages of Firmicutes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import os
#from rrnpp_detector.wrappers import orfipy
#from rrnpp_detector.integrate_annotations import orfipy_bed_to_dict
import orfipy_core
from rrnpp_detector.preprocessing import load_rbs_regex, load_translation_table
import re


def make_coordinates_to_annotated_peptides_dict(peptides, protein_dict):
    coordinates_to_annotated_peptides = dict()
    for p in peptides:
        genomic_key = protein_dict[p]['genomic_accession'] + '_' + protein_dict[p]['strand'] + '_' + str(protein_dict[p]['start']) + '_' + str(protein_dict[p]['end'])
        coordinates_to_annotated_peptides[genomic_key] = p
        protein_dict[p]['RBS_bin'] = protein_dict[p]['RBS_spacer'] = protein_dict[p]['RBS_motif'] = None
    return coordinates_to_annotated_peptides


def define_flanking_regions(anchor_proteins, protein_dict, position_dict, parameters):
    flanking_dict = dict()
    min_intergen_dist, max_intergen_dist = parameters['intergenic_distance_boundaries']
    max_small_orf_len = 3 * (parameters['propeptide_len_boundaries'][1] + 1)
    for anchor_protein in anchor_proteins:
        genomic_accession = protein_dict[anchor_protein]['genomic_accession']
        position = protein_dict[anchor_protein]['position']
        
        if not genomic_accession in flanking_dict:
            flanking_dict[genomic_accession] = dict()
        flanking_dict[genomic_accession][anchor_protein] = dict()
        
        # We want to search between receptor gene and flanking annotated gene
        # if flanking annotated gene too far away from receptor, then max_intergen_dist is used to set the boundary
        if (genomic_accession, position-1) in position_dict:
            annoprot_before = position_dict[(genomic_accession, position-1)]['protein_id']
            flanking_dict[genomic_accession][anchor_protein]['before_start'] = max((protein_dict[annoprot_before]['end'] - 20), (protein_dict[anchor_protein]['start'] - max_intergen_dist - max_small_orf_len))
        else:
            flanking_dict[genomic_accession][anchor_protein]['before_start'] = protein_dict[anchor_protein]['start'] - max_intergen_dist - max_small_orf_len
        flanking_dict[genomic_accession][anchor_protein]['before_end'] = protein_dict[anchor_protein]['start'] - (min_intergen_dist + 1)
        
        flanking_dict[genomic_accession][anchor_protein]['after_start'] = protein_dict[anchor_protein]['end'] + min_intergen_dist + 1
        if (genomic_accession, position+1) in position_dict:
            annoprot_after = position_dict[(genomic_accession, position+1)]['protein_id']
            flanking_dict[genomic_accession][anchor_protein]['after_end'] = min((protein_dict[annoprot_after]['start'] + 20), (protein_dict[anchor_protein]['end'] + max_intergen_dist + max_small_orf_len))
        else:
            flanking_dict[genomic_accession][anchor_protein]['after_end'] = protein_dict[anchor_protein]['end'] + max_intergen_dist + max_small_orf_len
    return flanking_dict


def search_rbs(dna_seq, rbs_dict):
    rbs_bin, rbs_motif, rbs_spacer = [0, '', None]
    for key in rbs_dict:
        match = rbs_dict[key]["regex"].search(dna_seq)
        if match:
            rbs_start = match.start()
            rbs_stop = rbs_start + rbs_dict[key]["length"]
            rbs_spacer = len(dna_seq) - rbs_stop
            rbs_motif = dna_seq[rbs_start:rbs_stop]
            rbs_bin = rbs_dict[key]['bin']
            break
    return rbs_bin, rbs_motif, rbs_spacer
                        
                        
def get_orfs(genomic_accession, strand, region_start, region_end, region_seq,
             min_orf_len, max_orf_len, start_codons, rbs_regex, bins_allowed, 
             orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
             pair_dict, cognate_dict, coord_to_anno_pep, protein_dict):
    orfs_coordinates = list()
    for start, stop, s, d in orfipy_core.orfs(region_seq, minlen = min_orf_len, maxlen = max_orf_len, starts = start_codons, strand = 'f'):
        orf_seq = region_seq[start:stop]
        
        orfs_coordinates.append([start, stop])
        n = len(orf_seq)
        # get nested orfs
        for i in range(3, n, 3):
            if orf_seq[i:3] in start_codons and (n-i) >= min_orf_len:
                orfs_coordinates.append([i+start, stop])
                
    for orf in orfs_coordinates:
        orf_id = 'smallORF' + str(orf_counter)
        start, stop = orf[0], orf[1]
        orf_seq = region_seq[start:stop]
        

        # search RBS
        upstream_start = start - 21
        upstream_stop = start
        if upstream_start < 0:
            upstream_start = 0
        upstream = region_seq[upstream_start:upstream_stop]
        rbs_bin, rbs_motif, rbs_spacer = search_rbs(upstream, rbs_regex)
        

        if strand == '-':
            real_start = region_end - stop + 1 - 3
            real_end = region_end - start
        else:
            real_start = start + region_start
            real_end = stop + region_start - 1 + 3 # (+3 for stop codon)
            
        genomic_key = genomic_accession + '_' + strand + '_' + str(real_start) + '_' + str(real_end)
        # if detected small orf already is annotated, don't print it in fasta
        # Store nonetheless the RBS values in protein dict
        if genomic_key in coord_to_anno_pep:
            continue
        
        # if orf has a relevant RBS, fill protein_dict, cognate_dict and pair_dict
        if rbs_bin in bins_allowed:
            protein_dict[orf_id] = dict()
            protein_dict[orf_id]['genomic_accession'] = genomic_accession
            protein_dict[orf_id]['strand'] = strand
            protein_dict[orf_id]['start'] = real_start
            protein_dict[orf_id]['end'] = real_end
            protein_dict[orf_id]['RBS_bin'] = rbs_bin
            protein_dict[orf_id]['RBS_spacer'] = rbs_spacer
            protein_dict[orf_id]['RBS_motif'] = rbs_motif
            pair_dict[pair_id] = [orf_id, anchor_protein]
            if not orf_id in cognate_dict:
                cognate_dict[orf_id] = dict()
                cognate_dict[orf_id]['cognate'] = list()
                cognate_dict[orf_id]['couple_id'] = list()
            cognate_dict[orf_id]['cognate'].append(anchor_protein)
            cognate_dict[orf_id]['couple_id'].append(pair_id)
            if not anchor_protein in cognate_dict:
                cognate_dict[anchor_protein] = dict()
                cognate_dict[anchor_protein]['cognate'] = list()
                cognate_dict[anchor_protein]['couple_id'] = list()
            cognate_dict[anchor_protein]['cognate'].append(orf_id)
            cognate_dict[anchor_protein]['couple_id'].append(pair_id)           
            outfile.write('>' + str(orf_id) + ' ' + genomic_accession + '_[' + str(real_start) + '-' + str(real_end) + '](' + strand + 
                          ')_Anchor=' + anchor_protein + '_RBS_bin=' + str(rbs_bin) + ';RBS_motif=' + rbs_motif + ';RBS_spacer=' + str(rbs_spacer) + '\n')
            outfile.write(orf_seq + '\n')
            non_annotated_peptides.append(orf_id)
            orf_counter += 1
            pair_id += 1
    return orf_counter, pair_id, non_annotated_peptides


def translate(input_fna, output_faa, writing_mode):
    translation_table = load_translation_table()
    with open(output_faa, mode=writing_mode) as outfile:
        with open(input_fna, mode='r') as infile:
            for line in infile:
                if line[0] == ">":
                    outfile.write(line)
                else:
                    protein = ''
                    cds = line.strip()
                    for i in range(0, len(cds), 3):
                        codon = cds[i:(i + 3)]
                        if codon in translation_table:
                            protein += translation_table[codon]
                        else:
                            protein += 'X'
                    for i in range(0, len(protein), 80):
                        outfile.write(protein[i:(i+80)] + '\n')
        infile.close()
    outfile.close()
        
# SPRAT stands for Small Peptides with RBS Annotation Tool      
def sprat(orf_counter, fna, out_faa, anchor_proteins, annotated_peptides, parameters, working_dir, protein_dict, position_dict, pair_dict, cognate_dict):
    if annotated_peptides:
        # create a dict to test if returned coordinates correspond to already existing peptides
        coord_to_anno_pep = make_coordinates_to_annotated_peptides_dict(annotated_peptides, protein_dict)
        pair_id = len(pair_dict)
    else:
        coord_to_anno_pep = dict()
        pair_id = 0
    non_annotated_peptides = list()
    # define flanking regions of candidate receptors
    flanking_dict = define_flanking_regions(anchor_proteins, protein_dict, cognate_dict, parameters)
    # sprat algo
    output_fna = os.path.join(working_dir, '1_possible_smORFs_flanking_candidate_receptors.fna')
    cdna = str.maketrans({'A':'T', 'C':'G', 'G':'C', 'T':'A'})
    min_orf_len = parameters['propeptide_len_boundaries'][0] * 3
    rbs_regex = load_rbs_regex()
    max_orf_len = parameters['propeptide_len_boundaries'][1] * 3
    is_bin_allowed = dict.fromkeys(parameters['allowed_rbs_bins'], True)
    with open(output_fna, mode='w') as outfile:
        dna_sequence, genomic_accession = ['', '']
        with open(fna, mode='r') as infile:
            for line in infile:
                if line[0] == ">":
                    if genomic_accession in flanking_dict:
                        anchor_proteins = flanking_dict[genomic_accession].keys()
                        for anchor_protein in anchor_proteins:
                            # region before anchor protein on positive strand
                            region_start = flanking_dict[genomic_accession][anchor_protein]['before_start']
                            region_end = flanking_dict[genomic_accession][anchor_protein]['before_end']
                            region_dna = dna_sequence[(region_start-1):region_end].upper()
                            orf_counter, pair_id, non_annotated_peptides = get_orfs(genomic_accession, '+', region_start, region_end, region_dna,
                                                                               min_orf_len, max_orf_len, parameters['orfipy_start_codons'], rbs_regex, is_bin_allowed, 
                                                                               orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
                                                                               pair_dict, cognate_dict, coord_to_anno_pep, protein_dict)
                            
                            # region before anchor protein on negative strand
                            region_dna = region_dna.translate(cdna)[::-1]
                            orf_counter, pair_id, non_annotated_peptides = get_orfs(genomic_accession, '-', region_start, region_end, region_dna,
                                                                               min_orf_len, max_orf_len, parameters['orfipy_start_codons'], rbs_regex, is_bin_allowed, 
                                                                               orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
                                                                               pair_dict, cognate_dict, coord_to_anno_pep, protein_dict)
                            
                            # region after anchor protein on positive strand
                            region_start = flanking_dict[genomic_accession][anchor_protein]['after_start']
                            region_end = flanking_dict[genomic_accession][anchor_protein]['after_end']
                            region_dna = dna_sequence[(region_start-1):region_end].upper()
                            orf_counter, pair_id, non_annotated_peptides = get_orfs(genomic_accession, '+', region_start, region_end, region_dna,
                                                                               min_orf_len, max_orf_len, parameters['orfipy_start_codons'], rbs_regex, is_bin_allowed, 
                                                                               orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
                                                                               pair_dict, cognate_dict, coord_to_anno_pep, protein_dict)
                            
                            # region after anchor protein on negative strand
                            region_dna = region_dna.translate(cdna)[::-1]
                            orf_counter, pair_id, non_annotated_peptides = get_orfs(genomic_accession, '-', region_start, region_end, region_dna,
                                                                               min_orf_len, max_orf_len, parameters['orfipy_start_codons'], rbs_regex, is_bin_allowed, 
                                                                               orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
                                                                               pair_dict, cognate_dict, coord_to_anno_pep, protein_dict)
                    genomic_accession = line[1:].strip().split(' ', 1)[0]
                    dna_sequence = ''
                else:
                    dna_sequence += line.strip()
        infile.close()
        # handle last regions
        if genomic_accession in flanking_dict:
            anchor_proteins = flanking_dict[genomic_accession].keys()
            for anchor_protein in anchor_proteins:
                # region before anchor protein on positive strand
                region_start = flanking_dict[genomic_accession][anchor_protein]['before_start']
                region_end = flanking_dict[genomic_accession][anchor_protein]['before_end']
                region_dna = dna_sequence[(region_start-1):region_end].upper()
                orf_counter, pair_id, non_annotated_peptides = get_orfs(genomic_accession, '+', region_start, region_end, region_dna,
                                                                   min_orf_len, max_orf_len, parameters['orfipy_start_codons'], rbs_regex, is_bin_allowed, 
                                                                   orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
                                                                   pair_dict, cognate_dict, coord_to_anno_pep, protein_dict)
                
                # region before anchor protein on negative strand
                region_dna = region_dna.translate(cdna)[::-1]
                orf_counter, pair_id, non_annotated_peptides = get_orfs(genomic_accession, '-', region_start, region_end, region_dna,
                                                                   min_orf_len, max_orf_len, parameters['orfipy_start_codons'], rbs_regex, is_bin_allowed, 
                                                                   orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
                                                                   pair_dict, cognate_dict, coord_to_anno_pep, protein_dict)
                
                # region after anchor protein on positive strand
                region_start = flanking_dict[genomic_accession][anchor_protein]['after_start']
                region_end = flanking_dict[genomic_accession][anchor_protein]['after_end']
                region_dna = dna_sequence[(region_start-1):region_end].upper()
                orf_counter, pair_id, non_annotated_peptides = get_orfs(genomic_accession, '+', region_start, region_end, region_dna,
                                                                   min_orf_len, max_orf_len, parameters['orfipy_start_codons'], rbs_regex, is_bin_allowed, 
                                                                   orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
                                                                   pair_dict, cognate_dict, coord_to_anno_pep, protein_dict)
                
                # region after anchor protein on negative strand
                region_dna = region_dna.translate(cdna)[::-1]
                orf_counter, pair_id, non_annotated_peptides = get_orfs(genomic_accession, '-', region_start, region_end, region_dna,
                                                                   min_orf_len, max_orf_len, parameters['orfipy_start_codons'], rbs_regex, is_bin_allowed, 
                                                                   orf_counter, pair_id, anchor_protein, outfile, non_annotated_peptides,
                                                                   pair_dict, cognate_dict, coord_to_anno_pep, protein_dict)
    outfile.close()
    writing_mode = 'w'
    if os.path.exists(out_faa):
        writing_mode = 'a'
    translate(output_fna, out_faa, writing_mode)
    return (orf_counter + 1), non_annotated_peptides
       
    
