#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-type quorum sensing systems in chromosomes, plasmids and phages of Firmicutes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import os
import pandas
import sys
import shutil
from math import ceil
from rrnpp_detector.filter import filter_by_length, filter_by_adjacency, subset_by_id, setdiff, get_only_partners
from rrnpp_detector.integrate_annotations import split_annotation_file, gff_to_dict, feature_tbl_to_dict, explicitize_duplicates
from rrnpp_detector.wrappers import hmmsearch, tprpred, signalp, predisi, blastp
from rrnpp_detector.search_small_orfs import sprat
from rrnpp_detector.make_summary import make_summary, merge_signalp_predisi, make_iterative_search_summary


def concatenate_files(input_files, output_file):
    with open(output_file, mode = 'wb') as outfile:
        for input_file in input_files:
            with open(input_file, mode = 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
            infile.close()
    outfile.close()

             
def append_files(old_file, new_file):
    with open(old_file, mode='a') as f1:
        with open(new_file, mode='r') as f2:
            for line in f2:
                f1.write(line)
        f2.close()
    f1.close()


def launcher(args, parameters, rrnpp_detector_dir, current_dir, out_dir, working_dir, fasta_dir, 
             fna, faa, annotation_file, annotation_format):
    
    ###############################################################################################
    # Some variables to initialize
    ###############################################################################################
    stop_codon_in_protein_seq = False
    search_small_orfs = False
    if fna:
        search_small_orfs = True # if fna is provided, a search for small orfs encoded in the vicinity of receptors will be performed
        if not faa:
            stop_codon_in_protein_seq = True # Prodigal adds the "*" symbol at the end of protein sequences
    initial_working_dir = working_dir
    if args.expand_to_homologs:
        query_receptors = list() 
        query_propeptides = list()
        query_propeptides_faa = os.path.join(initial_working_dir, "query_propeptides_for_iterative_search.faa")
        open(query_propeptides_faa, 'w').close()
    
    ###############################################################################################
    # Faa is filtered to retain only proteins of len compatible with that of receptors
    ###############################################################################################
    print('* Filtering faa by length compatible with that of receptors ...')
    same_len_as_receptors_faa = os.path.join(working_dir, '0_same_len_as_receptors.faa')
    same_len_as_receptors = filter_by_length(faa, same_len_as_receptors_faa, parameters['receptor_len_boundaries'], stop_codon_in_protein_seq)
    if not same_len_as_receptors:
        sys.exit('  End of execution: no proteins passed the filter')
    print('  %d compatible proteins' % len(same_len_as_receptors))
    
    ###############################################################################################
    # Faa is filtered to retain only annotated proteins of len compatible with that of propeptides
    # if args.fna: 
    #    it is not prohibitive to have no potential annotated propeptides at this stage 
    #    because non annotated small orfs encoded in the genomic vicinity of potential receptors 
    #    will be searched subsequently
    ###############################################################################################
    if not search_small_orfs:
        print('* Filtering faa by length compatible with that of propeptides ...')
    same_len_as_propeptides_faa = os.path.join(working_dir, '0_annotated_same_len_as_propeptides.faa')
    same_len_as_propeptides = filter_by_length(faa, same_len_as_propeptides_faa, parameters['propeptide_len_boundaries'], stop_codon_in_protein_seq)
    if not search_small_orfs:
        if not same_len_as_propeptides:
            sys.exit('  End of execution: no proteins passed the filter')
        print('  %d compatible proteins' % len(same_len_as_propeptides))
         
    ###############################################################################################
    # Determine whether contigs should be processed altogether or one at a time
    ###############################################################################################
    if args.chunk_size:
        # Split gff or feature table based on nucleotide sequence id (different genomes, chromosomes, plasmids etc...)
        print('* Split %s into chunks of %s genomes' % (annotation_format, args.chunk_size))
        list_contigs = split_annotation_file(annotation_file, annotation_format, working_dir, int(args.chunk_size))
        print('  %d contig(s) identified in the %s' % (len(list_contigs), annotation_format))
        max_i = len(list_contigs)
        increment = int(args.chunk_size)
         
    else:
        max_i = 1
        increment = 1
    
    smORF_counter = 0
    qss_counter = 0
    iteration = 0
    for i in range(0, max_i, increment):
        iteration += 1
        if args.chunk_size:
            contig = list_contigs[i]
            print('############################################################')
            print('* Processing chunk %d starting from genome %d/%d: %s' % (iteration, i+1, len(list_contigs), contig))
            working_dir = os.path.join(initial_working_dir, 'chunk_' + str(iteration))
            os.makedirs(working_dir)
        
        ###########################################################################################
        # Read coordinates of ORFs from annotation file and integrate them in python dictionaries   
        ###########################################################################################
        print('* Storing annotated ORFs into python dictionaries ...')
        if annotation_format == 'gff':
            if args.chunk_size:
                duplicates, position_dict, protein_dict = gff_to_dict(os.path.join(initial_working_dir, 'chunk_' + str(iteration) + '.gff'))
            else:
                duplicates, position_dict, protein_dict = gff_to_dict(annotation_file)
        elif annotation_format == 'feature_tbl':
            if args.chunk_size:
                duplicates, position_dict, protein_dict = feature_tbl_to_dict(os.path.join(initial_working_dir, 'chunk_' + str(iteration) + '.tsv'))
            else:
                duplicates, position_dict, protein_dict = feature_tbl_to_dict(annotation_file)
        
        ###########################################################################################
        # Several sanity checks:
        # 1. Ensure one protein = one unique identifier
        #    If duplicated ids in the original faa, a new fasta file is created with duplicates
        #    renamed according to a unique identifier: (<protein_id>__dupli_nb<n>)
        #    Duplicated identifiers can typically arise when the assembly has been downloaded 
        #    from the NCBI using RefSeq as the source database instead of GenBank
        # 2. Make sure protein id in faa = protein id in annotation file
        #    This test is made by extracting the id from the first header line ('^>') in faa
        #    and checking it can be use as a key of the protein dictionary
        # 3. If args.fna : make sure contig id in fna = contig id in annotation file
        ###########################################################################################
        if duplicates:
            same_len_as_receptors = explicitize_duplicates(same_len_as_receptors_faa, working_dir, duplicates)
            same_len_as_propeptides = explicitize_duplicates(same_len_as_propeptides_faa, working_dir, duplicates)       
        # Check faa and annotation intersection has been successful 
        # Namely, if the protein_id extracted from the seq header in faa is found in the dict of protein_ids extracted from the annotation file
        if i == 0:
            first_protein = list(protein_dict.keys())[0]
            try:
                subset_by_id(faa, os.path.join(working_dir, 'test_subsetting.faa'), [first_protein])
            except:
                sys.exit('error: intersection between faa and %s based on protein_id is not possible. '
                         'e.g. the protein_id \'%s\' extracted from the %s is not identified as a protein identifier in the faa' 
                         % (annotation_format, first_protein, annotation_format))
            if search_small_orfs:
                try:
                    subset_by_id(fna, os.path.join(working_dir, 'test_subsetting.fna'), [protein_dict[first_protein]['genomic_accession']])
                except:
                    print('warning: intersection between fna and %s based on genomic_accession is not possible. '
                          'e.g. the genomic_accession \'%s\' extracted from the %s is not identified as a sequence identifier in the fna\n'
                          'NB: the search for small ORFs will be disabled. Only annotated proteins will be used'
                          % (annotation_format, protein_dict[first_protein]['genomic_accession'], annotation_format))
                    search_small_orfs = False
        
        ###########################################################################################
        # If contigs processed one at a time: intersect list of encoded proteins with list of 
        #    annotated proteins of len compatible with that of receptors and propeptides
        # If small orf search possible: 
        #    it is not prohibitive to have no potential annotated propeptides at this stage
        ###########################################################################################
        if args.chunk_size:
            # Intersect list of all small proteins with list of proteins encoded by the current genetic element (same rationale for receptors)
            if same_len_as_propeptides:
                annotated_candidate_propeptides = list(set.intersection(set(list(protein_dict.keys())), set(same_len_as_propeptides)))
                if not search_small_orfs:
                    if not annotated_candidate_propeptides:
                        print('  Skip current chunk: no proteins of length compatible with that of propeptides')
                        continue
                propeptides_faa = os.path.join(working_dir, '0_same_len_as_propeptides_in_chunk.faa')
                subset_by_id(same_len_as_propeptides_faa, propeptides_faa, annotated_candidate_propeptides)
            candidate_receptors = list(set.intersection(set(list(protein_dict.keys())), set(same_len_as_receptors)))
            if not candidate_receptors:
                #print(same_len_as_receptors)
                print(protein_dict.keys())
                print('  Skip current chunk: no proteins of length compatible with that of receptors')
                continue
            receptors_faa = os.path.join(working_dir, '0_same_len_as_receptors_in_chunk.faa')
            subset_by_id(same_len_as_receptors_faa, receptors_faa, candidate_receptors)
        else:
            if same_len_as_propeptides:
                # annotated_candidate_propeptides = same_len_as_propeptides
                annotated_candidate_propeptides = [protein for protein in set(same_len_as_propeptides) if protein in protein_dict]
                if len(annotated_candidate_propeptides) < len(same_len_as_propeptides):
                    missing_proteins = setdiff(same_len_as_propeptides, annotated_candidate_propeptides)
                    print('warning: could not intersect annotation file and faa for the following proteins')
                    print(missing_proteins)
                    # that means some (not all) proteins in fasta cannot be found in the gff: these proteins are then removed from the search space
                    propeptides_faa = os.path.join(working_dir, '0_same_len_as_propeptides_in_protein_dict.faa')
                    subset_by_id(same_len_as_propeptides_faa, propeptides_faa, annotated_candidate_propeptides)
                    os.remove(same_len_as_propeptides_faa)
                    same_len_as_propeptides_faa = propeptides_faa
                else:
                    propeptides_faa = same_len_as_propeptides_faa
                if not search_small_orfs:
                    if not annotated_candidate_propeptides:
                        sys.exit('  End of execution: no proteins of length compatible with that of propeptides')        
            # candidate_receptors = same_len_as_receptors
            candidate_receptors = [protein for protein in set(same_len_as_receptors) if protein in protein_dict]
            if not candidate_receptors:
                sys.exit('  Skip current genetic element: no proteins of length compatible with that of receptors')
            if len(candidate_receptors) < len(same_len_as_receptors):
                # that means some (not all) proteins in fasta cannot be found in the gff: these proteins are then removed from the search space
                missing_proteins = setdiff(same_len_as_receptors, candidate_receptors)
                print('warning: could not intersect annotation file and faa for the following proteins')
                print(missing_proteins)
                receptors_faa = os.path.join(working_dir, '0_same_len_as_receptors_in_protein_dict.faa')
                subset_by_id(same_len_as_receptors_faa, receptors_faa, candidate_receptors)
                os.remove(same_len_as_receptors_faa)
                same_len_as_receptors_faa = receptors_faa
            else:
                receptors_faa = same_len_as_receptors_faa
            
        
        ###########################################################################################
        # If small orf search impossible: 
        #    reduce the search space to receptors adjacent to propeptides and vice-versa
        # if small orf search possible: 
        #    keep initial list of receptors 
        #    and update list of annotated propeptides based on adjacency with potential receptors
        ###########################################################################################
        if annotated_candidate_propeptides:
            if not search_small_orfs:
                print('* Selecting potential annotated propeptides and receptors based on coding sequence adjacency ...')
            annotated_candidate_propeptides, tmp_receptors, pair_dict, cognate_dict = filter_by_adjacency(annotated_candidate_propeptides, candidate_receptors, 0, protein_dict, position_dict)
            if not search_small_orfs:
                candidate_receptors = tmp_receptors
                if not annotated_candidate_propeptides or not candidate_receptors:
                    if args.chunk_size:
                        print('  Skip current chunk: no proteins passed the filter')
                        continue
                    else:
                        sys.exit('  End of execution: no proteins passed the filter')
                input_receptors_faa = receptors_faa
                receptors_faa = os.path.join(working_dir, '0_same_len_as_receptors_adj_to_annotated_same_len_as_propeptides.faa')
                subset_by_id(input_receptors_faa, receptors_faa, candidate_receptors)
            if annotated_candidate_propeptides:
                input_propeptides_faa = propeptides_faa
                propeptides_faa = os.path.join(working_dir, '0_annotated_same_len_as_propeptides_adj_to_same_len_as_receptors.faa')
                
                subset_by_id(input_propeptides_faa, propeptides_faa, annotated_candidate_propeptides)
            if not search_small_orfs:
                print('  %d potential propeptides' % len(annotated_candidate_propeptides))
                print('  %d potential receptors' % len(candidate_receptors))
        else:
            pair_dict = dict()
            cognate_dict = dict()
        
        
        ###########################################################################################
        # Filter receptors that do not contain TPRs
        # Distinguish between canonical RRNPP receptors and Rgg receptors based on TIGR01716 hmm
        ###########################################################################################
        print('* Selecting potential receptors with TPRs ...')
        print('  Running hmmsearch ...')
        tpr_df = hmmsearch(receptors_faa, 
                           os.path.join(rrnpp_detector_dir, 'data', 'hmm_libraries', 'TPRs.hmm'), 
                           os.path.join(working_dir, 'hmmsearch_TPRs_stdout.txt'),
                           os.path.join(working_dir, 'hmmsearch_TPRs_tblout.txt'),
                           os.path.join(working_dir, 'hmmsearch_TPRs_domtblout.txt'),
                           args.cpu, parameters['hmmsearch_max_evalue'])
        # Receptors from the Rgg-RopB family are characterized by the TIGR01716 hmm"
        candidate_rggs = tpr_df.loc[(tpr_df['hmm_name'] == 'TIGR01716')]['target_name'].tolist()
        candidate_receptors = setdiff(tpr_df['target_name'].tolist(), candidate_rggs)
        input_receptors_faa = receptors_faa
        receptors_faa = os.path.join(working_dir, '1_candidate_receptors_with_TPRs.faa')
        junk_faa = None
        if args.tprpred:
            if candidate_receptors or candidate_rggs:
                junk_faa = os.path.join(working_dir, '1_tprpred_targets.faa')
                subset_by_id(input_receptors_faa, receptors_faa, candidate_receptors + candidate_rggs, junk_faa)
            else:
                junk_faa = input_receptors_faa
        else:
            subset_by_id(input_receptors_faa, receptors_faa, candidate_receptors + candidate_rggs)
        if args.tprpred:
            print('  Running tprpred on proteins potentially missed by hmmsearch ...')
            tprpred_out = os.path.join(working_dir, 'tprpred_tblout.tsv')
            tprpred_df = tprpred(junk_faa, tprpred_out, rrnpp_detector_dir, current_dir, parameters['tprpred_min_probab'])
            candidate_receptors += tprpred_df['target_name'].tolist()
            if candidate_receptors:
                tpr_df = pandas.concat([tpr_df, tprpred_df])
                subset_by_id(input_receptors_faa, receptors_faa, candidate_receptors + candidate_rggs)
        if not candidate_receptors and not candidate_rggs:
            if args.chunk_size:
                print('  Skip current chunk: no potential receptors passed the filter')
                continue
            else:
                sys.exit('  End of execution: no proteins passed the filter')
        print('  %d potential receptors' % (len(candidate_receptors) + len(candidate_rggs)))
        search_space_receptors_faa = receptors_faa
        
        ###########################################################################################
        # Reduce the search space of annotated peptides to:
        #    - candidate propeptides adjacent to candidate receptors with TPRs
        #    - candidate SHPs adjacent to candidate Rgg receptors
        # if small orf search impossible:
        #    - candidate receptors without adjacent candidate propeptides/SHPs are filtered out
        ###########################################################################################
        if not search_small_orfs:
                print('* Selecting potential propeptides adjacent to receptors ...')
        tmp_receptors, annotated_candidate_propeptides = get_only_partners(candidate_receptors, cognate_dict)
        tmp_rggs, annotated_candidate_shps = get_only_partners(candidate_rggs, cognate_dict)
        if not search_small_orfs:
            if not annotated_candidate_propeptides and not annotated_candidate_shps:
                if args.chunk_size:
                    print('  Skip current chunk: no proteins passed the filter')
                    continue
                else:
                    sys.exit('  End of execution: no proteins passed the filter')
            candidate_receptors = tmp_receptors
            candidate_rggs = tmp_rggs
        input_propeptides_faa = propeptides_faa
        propeptides_faa = os.path.join(working_dir, '1_micropeptides_flanking_candidate_receptors_with_TPRs.faa')
        shps_faa = os.path.join(working_dir, '1_micropeptides_flanking_candidate_Rgg_receptors.faa')
        if annotated_candidate_propeptides:
            subset_by_id(input_propeptides_faa, propeptides_faa, annotated_candidate_propeptides)
        if annotated_candidate_shps:
            subset_by_id(input_propeptides_faa, shps_faa, annotated_candidate_shps)
        if not search_small_orfs:
            print('  %d potential propeptides' % (len(annotated_candidate_propeptides) + len(annotated_candidate_shps)))

        candidate_propeptides = annotated_candidate_propeptides
        candidate_shps = annotated_candidate_shps
        ###########################################################################################
        # If small orf search possible: SPRAT (Small Peptides with RBS Annotation Tool) 
        #     In addition to annotated peptides, unannotated protein-coding small ORFs 
        #     preceded by a SD-RBS are searched in the vicinity of candidate receptors
        #     (NB: if an unannotated ORF corresponds to an annotated ORF: the former is ignored)
        #     They are then translated into amino-acid sequences and added to the faa 
        #     already containing annotated proteins
        #     These unannotated peptides are also added to the pair_dict and cognate_dict dictionaries 
        #     At the end of the operation, candidate receptors with no adjacent candidate propeptides 
        #     are filtered out 
        ###########################################################################################
        if search_small_orfs:
            unannotated_candidate_propeptides = list()
            unannotated_candidate_shps = list()
            print('* Searching small ORFs with RBS in flanking regions of candidate receptors ...')
            if candidate_receptors:
                smORF_counter, unannotated_candidate_propeptides = sprat(smORF_counter, fna, propeptides_faa, candidate_receptors, annotated_candidate_propeptides, parameters, working_dir, protein_dict, position_dict, pair_dict, cognate_dict)
                candidate_propeptides += unannotated_candidate_propeptides                
            if candidate_rggs:
                smORF_counter, unannotated_candidate_shps = sprat(smORF_counter, fna, shps_faa, candidate_rggs, annotated_candidate_shps, parameters, working_dir, protein_dict, position_dict, pair_dict, cognate_dict)
                candidate_shps += unannotated_candidate_shps
            if not candidate_propeptides and not candidate_shps:
                if args.chunk:
                    print('  Skip current chunk: no potential propeptides identified')
                    continue
                else:
                    sys.exit('  End of execution: no potential propeptides identified')
            print('  %d non-annotated small ORFs with RBS motif adjacent to candidate receptors' % (len(unannotated_candidate_propeptides) + len(unannotated_candidate_shps)))
            print('  %d annotated small ORFs adjacent to candidate receptors' % (len(candidate_propeptides) + len(candidate_shps) - len(unannotated_candidate_propeptides) - len(unannotated_candidate_shps)))
            print('* Filtering receptors with no adjacent potential propeptides ...')
            candidate_propeptides, tmp_receptors = get_only_partners(candidate_propeptides, cognate_dict)
            candidate_receptors = list(set(candidate_receptors).intersection(set(tmp_receptors)))
            candidate_shps, tmp_rggs = get_only_partners(candidate_shps, cognate_dict)
            candidate_rggs = list(set(candidate_rggs).intersection(set(tmp_rggs)))
            input_receptors_faa = receptors_faa
            receptors_faa = os.path.join(working_dir, '1_candidate_receptors_with_TPRs_adjacent_to_micropeptides.faa')
            subset_by_id(input_receptors_faa, receptors_faa, candidate_receptors + candidate_rggs)
            print('  %d potential receptors' % (len(candidate_receptors) + len(candidate_rggs)))
        if candidate_shps and candidate_propeptides:
            search_space_propeptides_faa = os.path.join(working_dir, 'all_potential_propeptides.faa')
            concatenate_files([propeptides_faa, shps_faa], search_space_propeptides_faa)
        elif candidate_shps:
            search_space_propeptides_faa = shps_faa
        else:
            search_space_propeptides_faa = propeptides_faa
            
        ###########################################################################################
        # At this stage, all candidate receptors will eventually be printed 
        # either in the permissive or in the conservative output folders
        # The following lines are meant to characterize in more depth these putative receptors:
        #   - hmmsearch of DNA binding domains to identify transcription factors
        #   - blastp of reference receptors against candidate receptors
        ###########################################################################################
        print('* Characterizing potential receptors ...')
        print('  Transcription factor inference with hmmsearch ...')
        tf_df = hmmsearch(receptors_faa, 
                          os.path.join(rrnpp_detector_dir, 'data', 'hmm_libraries', 'HTH.hmm'), 
                          os.path.join(working_dir, 'hmmsearch_HTH_stdout.txt'),
                          os.path.join(working_dir, 'hmmsearch_HTH_tblout.txt'),
                          os.path.join(working_dir, 'hmmsearch_HTH_domtblout.txt'),
                          args.cpu, parameters['hmmsearch_max_evalue'])
        print('  %d potential receptors with a detected DNA binding domain' % len(tf_df))
        print('  Homology assessment to experimentally-validated receptors with Blastp ...')
        blast_df = blastp(os.path.join(rrnpp_detector_dir, 'data', 'fasta', 'ref_receptors.faa'), receptors_faa, working_dir, args.cpu, 
                          parameters['blastp_min_pident'], parameters['blastp_min_pcover'], parameters['blastp_max_evalue'])
        blast_df.rename(columns = {'query_name': 'homologous_ref_receptor'}, inplace = True)
        print('  %d potential receptors matched experimentally-validated receptors' % blast_df['homologous_ref_receptor'].count())
        
        ###########################################################################################
        # Candidate SHPs will be retained in the conservative output only if matching the SHP HMM
        # Candidate propeptides only if categorized by signalp or PrediSi as SP(Sec/SPI)
        ###########################################################################################
        # HMM search of SHP profile is searched against peptides adjacent to Rgg receptors
        # - Peptides in the neighborhood of Rgg receptors with no detected adjacent SHP
        #   are put in the pool of input peptides for signalP
        # - Peptides in the neighborhood of Rgg receptors with an adjacent SHP are discarded 
        #   (except of course SHP)
        ###########################################################################################
        print('* Selecting propeptides with a signal sequence or matching SHP profile ...')
        shps = list()
        permissive_shps = list()
        if candidate_shps:
            print('  Running hmmsearch to detect potential SHPs propeptides ...')
            tmp_shp_df = hmmsearch(shps_faa, 
                                   os.path.join(rrnpp_detector_dir, 'data', 'hmm_libraries', 'SHP.hmm'), 
                                   os.path.join(working_dir, 'hmmsearch_SHP_stdout.txt'),
                                   os.path.join(working_dir, 'hmmsearch_SHP_tblout.txt'),
                                   os.path.join(working_dir, 'hmmsearch_SHP_domtblout.txt'),
                                   args.cpu, parameters['hmmsearch_max_evalue'])
            shps = tmp_shp_df['target_name'].tolist()
            shps, rggs = get_only_partners(shps, cognate_dict)
            permissive_rggs = setdiff(candidate_rggs, rggs)
            permissive_rggs, permissive_shps = get_only_partners(permissive_rggs, cognate_dict)
            if permissive_shps:
                if candidate_propeptides:
                    s = set(candidate_propeptides)
                    permissive_shps = [protein for protein in permissive_shps if protein not in s]
                tmp_faa = os.path.join(working_dir, '1_Non-SHP_peptides_flanking_Rgg_receptors.faa')    
                subset_by_id(shps_faa, tmp_faa, permissive_shps)
            print('  %d detected SHPs' % len(shps))
            
        if permissive_shps:
            if not candidate_propeptides:
                propeptides_faa = tmp_faa
                candidate_propeptides = permissive_shps
            else:
                with open(tmp_faa, mode = 'r') as infile:
                    with open(propeptides_faa, mode = 'a') as outfile:
                        for line in infile:
                            outfile.write(line)
                    outfile.close()
                infile.close()
                candidate_propeptides += permissive_shps
                
        ###########################################################################################
        # SignalP is launched against peptides adjacent to receptors
        # - Peptides in the neighborhood of receptors with no detected adjacent propeptide
        #   are put in the pool of input peptides for Predisi
        # - Peptides in the neighborhood of receptors with an adjacent propeptide are discarded 
        #   (except of course peptides with a signal sequence)
        ###########################################################################################
        if candidate_propeptides:       
            print('  Running SignalP to predict secretion-mode of potential propeptides ...')
            print(propeptides_faa)
            signalp_df = signalp(propeptides_faa, working_dir, rrnpp_detector_dir, current_dir)
            propeptides = signalp_df.loc[signalp_df['Prediction'].isin(['SP(Sec/SPI)', 'LIPO(Sec/SPII)'])]['protein_id'].tolist()
            propeptides, receptors = get_only_partners(propeptides, cognate_dict)
            print('  %d detected propeptides with SignalP' % len(propeptides))
            
            #######################################################################################
            # predisi can be used to increase the sensitivity of the tool
            #######################################################################################
            if args.predisi and len(candidate_receptors) > len(receptors):
                print('  Running PrediSi on permissive candidates ...')
                permissive_receptors = setdiff(candidate_receptors, receptors)
                permissive_receptors, permissive_propeptides = get_only_partners(permissive_receptors, cognate_dict)
                predisi_input = os.path.join(working_dir, '1_predisi_targets.faa')
                subset_by_id(propeptides_faa, predisi_input, permissive_propeptides)
                predisi_df = predisi(predisi_input, working_dir, rrnpp_detector_dir, current_dir, parameters['predisi_min_likelihood'])
                if predisi_df is not None:
                    signalp_df = merge_signalp_predisi(signalp_df, predisi_df)
                    print('  %d detected propeptides with PrediSi' % len(predisi_df))
                    # propeptides += predisi_df['protein_id'].tolist()
                    # propeptides, receptors = get_only_partners(propeptides, cognate_dict)
                    # subset_by_id(input_propeptides_faa, propeptides_faa, shps + propeptides, junk_faa)
            
        ###########################################################################################
        # Append SHP results to Signalp_df
        ###########################################################################################
        if shps:
            # convert shp_df into same format as signalp_df and concatenate both dfs
            shp_df = pandas.DataFrame()
            shp_df['protein_id'] = tmp_shp_df['target_name']
            shp_df['Prediction'] = 'SHP_propeptide'
            shp_df['propeptide_score'] = tmp_shp_df['score'] * 0.01
            for c in signalp_df.columns[3:]:
                shp_df[c] = None
            if candidate_propeptides:
                signalp_df = pandas.concat([signalp_df, shp_df])
            else:
                signalp_df = shp_df
       
        ###########################################################################################
        # Merge receptors results (hmmsearch, blast) with propeptides results (signalp) 
        # based on cognate dictionary
        ###########################################################################################
        print('* Preparing output ...')
        tmp_df = make_summary(tpr_df, tf_df, blast_df, signalp_df, protein_dict, cognate_dict)
        tmp_df = tmp_df[tmp_df['intergenic_distance'].between(parameters['intergenic_distance_boundaries'][0], parameters['intergenic_distance_boundaries'][1])]
        
        ###########################################################################################
        # Strict output is composed of receptor-propeptide pairs with a propeptide categorized as 
        # - SHP according to hmmsearch or 
        # - SP(Sec/SPI) according to Signalp
        ###########################################################################################
        strict_df = tmp_df[(tmp_df['Prediction'].isin(parameters['strict_predictors'])) & (tmp_df['genomic_context'].isin(parameters['strict_contexts']))]
        strict_df = strict_df.assign(qss_id=range(qss_counter, qss_counter + len(strict_df)))
        strict_receptors = strict_df['receptor_id'].unique().tolist()
        strict_propeptides = strict_df['propeptide_id'].unique().tolist()
        qss_counter += len(strict_df)
        
        ###########################################################################################
        # Relaxed output is composed of receptor-propeptide pairs with neither a receptor or a 
        # propeptide present in the strict output and with a propeptide categorized as either 
        # - LIPO(Sec/SPII) according to Signalp or
        # - PrediSi_pos or PrediSi_neg_(but_high_score) according to PrediSi 
        ###########################################################################################
        tmp_df = tmp_df.loc[tmp_df['receptor_id'].isin(strict_receptors) == False]
        #relaxed_df = tmp_df.loc[(tmp_df['Prediction'].isin(relaxed_predictors)) | (tmp_df['SP(Sec/SPI)'] >= parameters['predisi_min_likelihood'])]
        relaxed_df = tmp_df.loc[tmp_df['Prediction'].isin(parameters['relaxed_predictors'])]
        relaxed_df = relaxed_df.assign(qss_id=range(qss_counter, qss_counter + len(relaxed_df)))
        relaxed_receptors = relaxed_df['receptor_id'].unique().tolist()
        relaxed_propeptides = relaxed_df['propeptide_id'].unique().tolist()
        qss_counter += len(relaxed_df)
        
        ###########################################################################################
        # Loose output is composed of receptor-propeptide pairs with neither a receptor or a 
        # propeptide present in the strict output or in the relaxed output
        # - For each receptor, only one peptide is outputted: the one with the highest RBS bin
        # - the only exception is when another neighboring peptide is predicted as TAT(Tat/SPI)
        ###########################################################################################
        tmp_df = tmp_df.loc[tmp_df['receptor_id'].isin(relaxed_receptors) == False]
        # tat_df = tmp_df.loc[tmp_df['Prediction'] == 'TAT(Tat/SPI)']
        # best_rbs_df = tmp_df.loc[tmp_df['Prediction'] != 'TAT(Tat/SPI)'].sort_values(by='p_RBS_bin', ascending=False, kind='mergesort').drop_duplicates(subset=['receptor_id'])
        loose_df_a = tmp_df.loc[tmp_df.propeptide_id.str[0:8] != "smallORF"]
        loose_df_b = tmp_df.loc[tmp_df['p_RBS_bin'].isin(parameters['most_used_rbs_bins'])].sort_values(by='p_RBS_bin', ascending=False, kind='mergesort').drop_duplicates(subset=['receptor_id'])
        loose_df = pandas.concat([loose_df_a, loose_df_b])
        # loose_df = pandas.concat([tat_df, best_rbs_df])
        loose_df = loose_df.assign(qss_id=range(qss_counter, qss_counter + len(loose_df)))
        loose_receptors = loose_df['receptor_id'].unique().tolist()
        loose_propeptides = loose_df['propeptide_id'].unique().tolist()
        qss_counter += len(loose_df)
        
        if contig:
            suffix_faa = "__chunk_" + str(iteration)
        else:
            suffix_faa = ""
        if strict_receptors:
            if not os.path.exists(os.path.join(out_dir, '1_summary_strict_output.tsv')):
                writing_mode, write_header = 'w', True
            else:
                writing_mode, write_header = 'a', False
            strict_df.to_csv(os.path.join(out_dir, '1_summary_strict_output.tsv'), sep='\t', header=write_header, index=False, mode=writing_mode)
            subset_by_id(search_space_receptors_faa, os.path.join(fasta_dir, "1_receptors__strict_detection" + suffix_faa + ".faa"), strict_receptors)
            subset_by_id(search_space_propeptides_faa, os.path.join(fasta_dir, "1_propeptides__strict_detection" + suffix_faa + ".faa"), strict_propeptides)
            if args.expand_to_homologs:
                query_receptors += strict_receptors
                query_propeptides += query_propeptides
                append_files(query_propeptides_faa, os.path.join(fasta_dir, "1_propeptides__strict_detection" + suffix_faa + ".faa"))
            print('* %d detected QSSs with detection_strictness = strict' % len(strict_df))
            
        if relaxed_receptors:
            if not os.path.exists(os.path.join(out_dir, '2_summary_relaxed_output.tsv')):
                writing_mode, write_header = 'w', True
            else:
                writing_mode, write_header = 'a', False
            relaxed_df.to_csv(os.path.join(out_dir, '2_summary_relaxed_output.tsv'), sep='\t', header=write_header, index=False, mode=writing_mode)
            subset_by_id(search_space_receptors_faa, os.path.join(fasta_dir, "2_receptors__relaxed_detection" + suffix_faa + ".faa"), relaxed_receptors)
            subset_by_id(search_space_propeptides_faa, os.path.join(fasta_dir, "2_propeptides__relaxed_detection" + suffix_faa + ".faa"), relaxed_propeptides)
            print('* %d detected QSSs with detection_strictness = relaxed' % len(relaxed_df))
                
        if loose_receptors:
            if not os.path.exists(os.path.join(out_dir, '3_summary_loose_output.tsv')):
                writing_mode, write_header = 'w', True
            else:
                writing_mode, write_header = 'a', False
            loose_df.to_csv(os.path.join(out_dir, '3_summary_loose_output.tsv'), sep='\t', header=write_header, index=False, mode=writing_mode) 
            print('* %d hypothetical QSSs with detection_strictness = loose' % len(loose_df))
            subset_by_id(search_space_receptors_faa, os.path.join(fasta_dir, "3_hypothetical_receptors__loose_detection" + suffix_faa + ".faa"), loose_receptors)
            subset_by_id(search_space_propeptides_faa, os.path.join(fasta_dir, "3_hypothetical_propeptides__loose_detection" + suffix_faa + ".faa"), loose_propeptides)
            
        if strict_receptors or relaxed_receptors or loose_receptors:
            strict_df.insert(0,  'detection_strictness', 'strict')
            relaxed_df.insert(0, 'detection_strictness', 'relaxed')
            loose_df.insert(0,   'detection_strictness', 'loose')
            extensive_df = pandas.concat([strict_df, relaxed_df, loose_df])
            if not os.path.exists(os.path.join(out_dir, '4_summary_all_outputs.tsv')):
                writing_mode, write_header = 'w', True
            else:
                writing_mode, write_header = 'a', False
            extensive_df.to_csv(os.path.join(out_dir, '4_summary_all_outputs.tsv'), sep='\t', header=write_header, index=False, mode=writing_mode) 
                  
    ###############################################################################################
    # Iterative blastp search
    ###############################################################################################
    if args.expand_to_homologs and query_receptors:
        print('############################################################')
        print("* Expansion to undetected homologs")
        if args.chunk_size:
            sys.exit("  Expansion to homologs is not yet available when target genomes are not processed altogether")
        query_receptors_faa = os.path.join(initial_working_dir, "query_receptors_for_iterative_search.faa")
        target_receptors_faa = os.path.join(initial_working_dir, "target_blastdb_of_receptors_for_iterative_search.faa")
        subset_by_id(same_len_as_receptors_faa, query_receptors_faa, query_receptors, target_receptors_faa)
        
        r_blast_df = blastp(query_receptors_faa, target_receptors_faa, working_dir, args.cpu, parameters['blastp_min_pident_iterative_search_receptors'], 
                            parameters['blastp_min_pcover_iterative_search_receptors'], parameters['blastp_max_evalue_iterative_search_receptors'])
        homologous_receptors = r_blast_df['target_name'].tolist()
        if not homologous_receptors:
            sys.exit('  End of execution: No homologous receptors detected ...')
        print("  %d receptors homologous to receptors detected with detection stritness = strict" % len(homologous_receptors))
        
        
        parameters['allowed_rbs_bins'] = [i for i in range(27, -1, -1)]
        tmp_receptors, annotated_candidate_propeptides = get_only_partners(homologous_receptors, same_len_as_propeptides)

        if annotated_candidate_propeptides:
            annotated_target_propeptides_faa = os.path.join(initial_working_dir, "annotated_propeptides_adj_to_homologous_receptors.faa")
            subset_by_id(same_len_as_propeptides_faa, annotated_target_propeptides_faa, annotated_candidate_propeptides)
        unannotated_candidate_propeptides = list()
        if search_small_orfs:
            unannotated_target_propeptides_faa = os.path.join(initial_working_dir, "unannotated_propeptides_adj_to_homologous_receptors.faa")
            print('* Searching small ORFs with RBS in flanking regions of homologous receptors ...')
            smORF_counter, unannotated_candidate_propeptides = sprat(smORF_counter, fna, unannotated_target_propeptides_faa, homologous_receptors, annotated_candidate_propeptides, parameters, working_dir, protein_dict, position_dict, pair_dict, cognate_dict)
        
        if annotated_candidate_propeptides and unannotated_candidate_propeptides:
            propeptides_faa = os.path.join(initial_working_dir, "all_propeptides_adj_to_homologous_receptors.faa")
            concatenate_files([annotated_target_propeptides_faa, unannotated_target_propeptides_faa], propeptides_faa)
        elif annotated_candidate_propeptides:
            propeptides_faa = annotated_target_propeptides_faa
        elif unannotated_candidate_propeptides:
            propeptides_faa = unannotated_target_propeptides_faa
        else:
            sys.exit('  End of execution: No adjacent peptides identified ...')    
        
        p_blast_df = blastp(query_propeptides_faa, propeptides_faa, working_dir, args.cpu, parameters['blastp_min_pident_iterative_search_propeptides'], 
                            parameters['blastp_min_pcover_iterative_search_propeptides'], parameters['blastp_max_evalue_iterative_search_propeptides'])
        homologous_propeptides = p_blast_df['target_name'].tolist()
        if not homologous_propeptides:
            sys.exit('  End of execution: No homologous propeptides detected ...')
        print("  %d adjacent propeptides homologous to propeptides detected with detection stritness = strict" % len(homologous_propeptides))
        
        expansion_dir = os.path.join(out_dir, 'expansion_to_homologs')
        os.makedirs(expansion_dir)
        expansion_df = make_iterative_search_summary(r_blast_df, p_blast_df, protein_dict, cognate_dict)    
        expansion_df.to_csv(os.path.join(expansion_dir, 'undetected_homologous_systems.tsv'), sep='\t', header=True, index=False, mode='w')
        
        subset_by_id(propeptides_faa, os.path.join(expansion_dir, "homologuous_propeptides_adjacent_to_homologous_receptors.faa"), expansion_df['homologous_propeptide'].unique().tolist())
        subset_by_id(target_receptors_faa, os.path.join(expansion_dir, "homologuous_receptors_adjacent_to_homologous_propeptides.faa"), expansion_df['homologous_receptor'].unique().tolist())

    if not args.keep_working_dir:
        shutil.rmtree(initial_working_dir)
    print('Completed!')
