#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in gram-positive bacteria and bacteriophages genomes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import argparse
import os
import pandas
import shutil
import sys


def __main__():
    # Read arguments
    parser = argparse.ArgumentParser(description='RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in gram-positive bacteria and bacteriophages genomes')
    parser.add_argument('-o', dest='out_dir', type=str, default='./', help='path to the output directory (default is current directory)')
    parser.add_argument('--fna', dest='fna', type=str, help='path to the fasta of the target genome(s) (will run Prodigal to detect ORFs)')
    parser.add_argument('--faa', dest='faa', type=str, help='path to the fasta of the protein sequences of the target genome(s) (requires additional --gff or --ft option)')
    parser.add_argument('--gff', dest='gff', type=str, help='path to the annotations of the target genome(s) in gff')
    parser.add_argument('--ft', dest='feature_tbl', type=str, help='path to the annotations of the target genome(s) in the NCBI_assembly feature_table format')
    parser.add_argument('--cpu', dest='cpu', type=str, default='1', help='number of cpu to use (default is 1)')
    parser.add_argument('--preserve_ram', action='store_true', help='minimize RAM usage at the expense of speed (will process target genomes one by one instead of all together)')
    args = parser.parse_args()    
    rrnpp_detector_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    current_dir = os.path.abspath(os.getcwd())      
    
    # Load built-in libraries after argparse to ensure fast printing of the script usage (-h)
    from python import check_args, handle_annotations, handle_fasta, hmm, blast, signalp, make_summary
    
    # Load built-in parameters
    # I might consider enabling the user to change default parameters in the future
    propeptide_length_boundaries = [10, 70]
    receptor_length_boundaries = [250, 500]
    hmmsearch_max_evalue = 1E-5
    max_intergenic_distance = 400 # nt
    blastp_min_pident = 20
    blastp_min_pcover = 60
    blastp_max_evalue = 1E-5
    
    # Check args and files provided by the user are OK
    out_dir, fasta, fasta_type = check_args.check_args(args)
    
    # Define working directory (will eventually be deleted)
    working_dir = os.path.join(out_dir, 'working_dir')
    os.makedirs(working_dir)
    # Define fasta output directoriy
    fasta_out_dir = os.path.join(out_dir, 'fastas')
    os.makedirs(fasta_out_dir)

    # If fna provided, run prodigal to predict the ORFs (will output a faa + a gff)
    if fasta_type == 'fna':
        annotation_format = 'gff'
        print('-------------------------------------')
        print('* Running Prodigal to predict ORFs in the target genome ...')
        annotation_file, faa = handle_annotations.prodigal(fasta, working_dir)
    # Else use the annotation file and the faa provided by the user
    elif fasta_type == 'faa':
        faa = args.faa
        # Check annotation file complies with requirements of its format
        if args.feature_tbl:
            annotation_format = 'feature_tbl'
            annotation_file = args.feature_tbl
        elif args.gff:
            annotation_format = 'gff'
            annotation_file = args.gff
        handle_annotations.test_annotation_file(annotation_file, annotation_format, working_dir)

              
    # Step 0: reduce the search space to proteins with a length compatible with reference propeptides and receptors, respectively
    print('-------------------------------------')
    print('* Filtering faa by lengths compatible with propeptides (%d-%daa) and receptors (%d-%daa) ...' % 
          (propeptide_length_boundaries[0], propeptide_length_boundaries[1], receptor_length_boundaries[0], receptor_length_boundaries[1]))
    if fasta_type == 'fna':
        # Add +1 to sequence length boundaries because prodigal uses the '*' symbol for stop codon
        propeptide_length_boundaries[0] += 1
        propeptide_length_boundaries[1] += 1
        receptor_length_boundaries[0] += 1
        receptor_length_boundaries[1] += 1
    all_small_proteins = handle_fasta.filter_by_length(faa, os.path.join(working_dir, 'step0_propeptides.faa'), propeptide_length_boundaries[0], propeptide_length_boundaries[1])
    all_medium_proteins = handle_fasta.filter_by_length(faa, os.path.join(working_dir, 'step0_receptors.faa'), receptor_length_boundaries[0], receptor_length_boundaries[1])
    if not (all_small_proteins or all_medium_proteins):
        sys.exit('  End of execution: no proteins passed the filter')
    print('  Filter done!')
    
    if args.preserve_ram:
        # Split gff or feature table based on nucleotide sequence id (different genomes, chromosomes, plasmids etc...)
        print('-------------------------------------')
        print('* Split %s into as many files as there are target genetic elements (genomes, chromosomes, plasmids ...)' % annotation_format)
        list_seq_ids = handle_annotations.split_annotation_file(annotation_file, annotation_format, working_dir)
        print(' %d genetic element(s) identified in the %s' % (len(list_seq_ids), annotation_format))
        max_iterations = len(list_seq_ids)  
    else:
        max_iterations = 1
 
        
    # Number of genetic element in which at least one exploratory QSS has been found
    n_successful_genetic_element = 0
    
    # Do the next steps of the algorithm for each genetic element (except for fast mode (all genetic element processed at once)) 
    for i in range(0, max_iterations):
        
        if args.preserve_ram:
            seq_id = list_seq_ids[i]
            suffix_faa = '_' + seq_id
            print('#####################################')
            print('* Processing genetic element %d/%d: %s' % (i+1, max_iterations, seq_id))
        else:
            seq_id = suffix_faa = ''

            
        # Integrate genomic info of ORFs into python dictionaries (protein_dict and position_dict)
        print('-------------------------------------')
        print('* Storing annotated ORFs into python dictionaries ...')
        if annotation_format == 'gff':
            if args.preserve_ram:
                duplicates, position_dict, protein_dict = handle_annotations.gff_to_dict(os.path.join(working_dir, seq_id + '.gff'))
            else:
                duplicates, position_dict, protein_dict = handle_annotations.gff_to_dict(annotation_file)
            print('  Gff to dictionaries done!')
        elif annotation_format == 'feature_tbl':
            if args.preserve_ram:
                duplicates, position_dict, protein_dict = handle_annotations.feature_tbl_to_dict(os.path.join(working_dir, seq_id + '.txt'))
            else:
                duplicates, position_dict, protein_dict = handle_annotations.feature_tbl_to_dict(annotation_file)
            print('  Feature table to dictionaries done!')

        # Copy fasta but duplicate sequence that appears multiple times in the annotation file (multiple copies or identical genes across genomes) and assign them a specific header (<protein_id>__nb_dupli<n>)
        # It might typically be the case if the assembly has been downloaded from the NCBI using RefSeq as the source database instead of GenBank
        if len(duplicates) > 0:
            all_small_proteins = handle_fasta.explicitize_duplicates(os.path.join(working_dir, 'step0_propeptides.faa'), working_dir, duplicates)
            all_medium_proteins = handle_fasta.explicitize_duplicates(os.path.join(working_dir, 'step0_receptors.faa'), working_dir, duplicates)
            
        # Check faa and gff intersection has been successful 
        # Namely, if the protein_id extracted from the seq header in faa is found in the dict of protein_ids extracted from the annotation file
        if i == 0:
            first_protein = list(protein_dict.keys())[0]
            try:
                handle_fasta.subset_by_id(faa, os.path.join(working_dir, 'test_subsetting.faa'), [first_protein])
            except:
                sys.exit('error: intersection between faa and %s based on protein_id is not possible. e.g. the protein_id \'%s\' extracted from faa is not identified as a protein_id in the %d' % (first_protein, annotation_format))
        
        if args.preserve_ram:
            # Intersect list of all small proteins with list of proteins encoded by the current genetic element (same rationale for receptors)
            propeptides_ids = list(set.intersection(set(list(protein_dict.keys())), set(all_small_proteins)))
            receptors_ids = list(set.intersection(set(list(protein_dict.keys())), set(all_medium_proteins)))   
        else:
            propeptides_ids = all_small_proteins
            receptors_ids = all_medium_proteins
     
        if len(propeptides_ids) == 0 or len(receptors_ids) == 0:
            if args.preserve_ram:
                print('  Skip current genetic element: no proteins of length compatible with propeptides or receptors')
                continue  
            else:
                sys.exit('  End of execution: no proteins of length compatible with propeptides or receptors')        

        # Step 1: reduce the search space to receptors adjacent to propeptides and vice-versa
        print('-------------------------------------')
        print('* Filtering potential propeptides and receptors based on coding sequence adjacency ...')
        propeptides_ids, receptors_ids = handle_fasta.filter_by_adjacency(propeptides_ids, receptors_ids, protein_dict, position_dict)
        if not (propeptides_ids or receptors_ids):
            if args.preserve_ram:
                print('  Skip current genetic element: no proteins passed the filter')
                continue
            else:
                sys.exit('  End of execution: no proteins passed the filter')
        handle_fasta.subset_by_id(os.path.join(working_dir, 'step0_propeptides.faa'), os.path.join(working_dir, 'step1_propeptides.faa'), propeptides_ids)
        handle_fasta.subset_by_id(os.path.join(working_dir, 'step0_receptors.faa'), os.path.join(working_dir, 'step1_receptors.faa'), receptors_ids)
        print('  Filter done!')
        print('  %d potential propeptides at this stage of the analysis' % len(propeptides_ids))
        print('  %d potential receptors at this stage of the analysis' % len(receptors_ids))
    
        # Step 2: filter receptors that do not contain TPRs and retain only their adjacent propeptides
        print('-------------------------------------')
        print('* Filtering potential receptors that do not match HMMs of TPRs ...')
        print('  Running hmmsearch ...')
        receptors_ids, tpr_hmm_df = hmm.hmmsearch(os.path.join(working_dir, 'step1_receptors.faa'), os.path.join(rrnpp_detector_dir, 'HMMs', 'TPRs.hmm'),
                                                  os.path.join(working_dir, 'hmmsearch_TPRs_stdout.txt'), os.path.join(working_dir, 'hmmsearch_TPRs_tblout.txt'), args.cpu, hmmsearch_max_evalue)
        if not receptors_ids:
            if args.preserve_ram:
                print('  Skip current genetic element: no potential receptors passed the filter')
                continue
            else:
                sys.exit('  End of execution: no proteins passed the filter')
        print('  Filter done!')
        print('  %d potential receptors at this stage of the analysis' % len(receptors_ids))
        handle_fasta.subset_by_id(os.path.join(working_dir, 'step1_receptors.faa'), os.path.join(working_dir, 'step2_receptors.faa'), receptors_ids)
        print('-------------------------------------')
        print('* Filtering potential propeptides not adjacent to TPR-containing receptors ...')
        propeptides_ids, receptors_ids = handle_fasta.filter_by_adjacency(propeptides_ids, receptors_ids, protein_dict, position_dict)
        if not propeptides_ids:
            if args.preserve_ram:
                print('  Skip current genetic element: no potential propeptides passed the filter')
                continue
            else:
                sys.exit('  End of execution: no potential propeptides passed the filter')
        handle_fasta.subset_by_id(os.path.join(working_dir, 'step1_propeptides.faa'), os.path.join(working_dir, 'step2_propeptides.faa'), propeptides_ids)
        print('  %d potential propeptides at this stage of the analysis' % len(propeptides_ids))
        print('  %d potential receptors at this stage of the analysis' % len(receptors_ids))
        
        # Characterization steps
        # Identify receptors that might be transcription factors (containing a DNA binding domain)
        print('-------------------------------------')
        print('* Characterizating potential receptors ...')
        print('  Detecting potential receptors that match HMMs of DNA binding domains')
        print('  Running hmmsearch ...')
        tf, tf_hmm_df = hmm.hmmsearch(os.path.join(working_dir, 'step2_receptors.faa'), os.path.join(rrnpp_detector_dir, 'HMMs', 'DNA_binding_motifs.hmm'),
                                      os.path.join(working_dir, 'hmmsearch_DNA-binding_stdout.txt'), os.path.join(working_dir, 'hmmsearch_DNA-binding_tblout.txt'), args.cpu, hmmsearch_max_evalue)
        print('  Detection done!')
        print('  %d potential receptors are predicted transcription factors' % len(tf))
        
        # Identify receptors that are homologous to experimentally-validated receptors
        print('  Blasting reference receptors against potential receptors ...')
        blast_df = blast.blastp(os.path.join(working_dir, 'step2_receptors.faa'), working_dir, args.cpu, 
                                rrnpp_detector_dir, blastp_min_pident, blastp_min_pcover, blastp_max_evalue)
        print('  Blast done!')
        print('  %d potential receptors are homologous to experimentally-validated receptors' % blast_df['homologous_ref_receptor'].count())
        
        # Run signalP to predict the secretion mode of adjacent propeptides
        print('-------------------------------------')
        print('* Running SignalP to predict secretion-mode of potential propeptides ...')
        signalp_df = signalp.signalp(os.path.join(working_dir, 'step2_propeptides.faa'), working_dir, rrnpp_detector_dir, current_dir)
        print('  SignalP done!')
        print('-------------------------------------')
        
        # Merge receptors results (hmmsearch, blast) with propeptides results (signalp) based on coding sequence adjacency
        print('* Summarizing intermediate results (exploratory output) ...')
        exploratory_summary_df = make_summary.make_summary(tpr_hmm_df, tf_hmm_df, blast_df, signalp_df, protein_dict, position_dict)
        print('  %d potential quorum sensing systems at this stage of the analysis' % len(exploratory_summary_df))
        
        # Pandas.write_csv parameters for outputting the summaries
        if n_successful_genetic_element == 0:
            writing_mode = 'w'
            write_header = True   
        elif n_successful_genetic_element == 1:
            # Append mode
            writing_mode = 'a' 
            write_header = False
        n_successful_genetic_element += 1
        
        # Step 3: filter propeptides not predicted to be secreted via the SP(Sec/SPI) mode and Write output files
        print('-------------------------------------')
        print('* Filtering potential propeptides not predicted to be secreted via the SP(Sec/SPI) mode ...')
        conservative_summary_df = exploratory_summary_df.loc[(exploratory_summary_df['Prediction'] == 'SP(Sec/SPI)') & (exploratory_summary_df['intergenic_distance'] <= max_intergenic_distance)]
        conservative_propeptides_ids = conservative_summary_df['propeptide_id'].unique().tolist()
        conservative_receptors_ids = conservative_summary_df['receptor_id'].unique().tolist()
        if conservative_propeptides_ids:
            print('  %d potential propeptides at this stage of the analysis' % len(conservative_propeptides_ids))
            print('  %d potential receptors at this stage of the analysis' % len(conservative_receptors_ids))            
            print('-------------------------------------')
            print('* Summarizing final results (conservative output) ...')
            print('  %d potential quorum sensing systems detected' % len(conservative_summary_df))
            print('-------------------------------------')
            print('* Writing output files ...')
            if len(exploratory_summary_df) > len(conservative_summary_df):
                # Filter all qss in exploratory_summary_df that are present in conservative_summary_df
                only_exploratory_qss_ids = make_summary.setdiff(exploratory_summary_df['qss_id'].tolist(), conservative_summary_df['qss_id'].tolist())
                exploratory_summary_df = exploratory_summary_df[exploratory_summary_df['qss_id'].isin(only_exploratory_qss_ids)]
                exploratory_summary_df.drop('qss_id', axis=1).to_csv(os.path.join(out_dir, 'exploratory_summary.tsv'), sep='\t', header=write_header, index=False, mode=writing_mode)
                # The exploratory proteins are all the proteins that have been filtered out at step 3 
                exploratory_propeptides_ids = make_summary.setdiff(propeptides_ids, conservative_propeptides_ids)
                exploratory_receptors_ids = make_summary.setdiff(receptors_ids, conservative_receptors_ids)
                # Write the exploratory output faa
                if len(exploratory_propeptides_ids) > 0:
                    handle_fasta.subset_by_id(os.path.join(working_dir, 'step2_propeptides.faa'), os.path.join(fasta_out_dir, 'exploratory_candidate_propeptides' + suffix_faa + '.faa'), exploratory_propeptides_ids)
                if len(exploratory_receptors_ids) > 0:
                    handle_fasta.subset_by_id(os.path.join(working_dir, 'step2_receptors.faa'), os.path.join(fasta_out_dir, 'exploratory_candidate_receptors' + suffix_faa + '.faa'), exploratory_receptors_ids)
            # Write the conservative output files
            handle_fasta.subset_by_id(os.path.join(working_dir, 'step2_propeptides.faa'), os.path.join(fasta_out_dir, 'conservative_candidate_propeptides' + suffix_faa + '.faa'), conservative_propeptides_ids)
            handle_fasta.subset_by_id(os.path.join(working_dir, 'step2_receptors.faa'), os.path.join(fasta_out_dir, 'conservative_candidate_receptors' + suffix_faa + '.faa'), conservative_receptors_ids)
            conservative_summary_df.drop('qss_id', axis=1).to_csv(os.path.join(out_dir, 'conservative_summary.tsv'), sep='\t', header=write_header, index=False, mode=writing_mode)
            print('  Output files written!')
        else:
            print('  No potential propeptides passed the filter')
            print('-------------------------------------')
            print('* Writing exploratory results before analyzing next genetic element ...')
            shutil.move(os.path.join(working_dir, 'step2_propeptides.faa'), os.path.join(fasta_out_dir, 'exploratory_candidate_propeptides' + suffix_faa + '.faa'))
            shutil.move(os.path.join(working_dir, 'step2_receptors.faa'), os.path.join(fasta_out_dir, 'exploratory_candidate_receptors' + suffix_faa + '.faa'))
            exploratory_summary_df.drop('qss_id', axis=1).to_csv(os.path.join(out_dir, 'exploratory_summary.tsv'), sep='\t', header=write_header, index=False, mode=writing_mode)
            print('  Output files written!')
        
        print('-------------------------------------')
        print('Completed!')
        
    #shutil.rmtree(working_dir)
    
if __name__=='__main__': __main__()

