
import argparse
import copy
import datetime
import os
import pickle
import sys
import time

from colorama import Fore, Style    
from utilities import *

def main(args):
    sys.stdout.write(Fore.GREEN + '[MAIN PROGRAM]\n' + Style.RESET_ALL)
    start_time = time.time()
    sys.stdout.write(Fore.GREEN + '[PROCESS]' + Style.RESET_ALL + ' argument validation\n')
    input_file, args = process_and_validate_argument(args)
    
    if args.reference:
        sys.stdout.write(Fore.MAGENTA + '[REFERENCE]\n')
        if 'ref_read_pickle' not in input_file:
            sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' object construction\n')
            ref_seq_dict = construct_sequence(input_file['ref_fasta'])
            #if 'ref_xprs_file' in input_file:
            sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' import expression value\n')
            read_expression(input_file, ref_seq_dict, 'ref')
            sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' calculation of sequence ambiguity\n')
            ref_self_blastn = filter_blastn(input_file['ref_blastn'])
            ref_self_match_dict = intersect_match(ref_self_blastn, ref_seq_dict, copy.deepcopy(ref_seq_dict))
            sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' construction of connected-component based graph\n')
            ref_uf, ref_component_dict = construct_graph(ref_seq_dict, ref_self_match_dict)
            sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' cluster expression\n')
            ref_cluster_dict = cluster_xprs(input_file, ref_seq_dict, ref_component_dict, 'ref')
            if 'ref_transrate' in input_file:
                sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' import transrate score\n')
                read_transrate(input_file, ref_seq_dict, 'ref')
            if 'ref_gtf' in input_file:
                sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' gene-isoform matching\n')
                gene_isoform_analysis(input_file['ref_gtf'], ref_seq_dict)   
            if 'ref_write_pickle' in input_file:   
                with open(input_file['ref_write_pickle'], 'wb') as ref_write_pickle:
                    pickle.dump((ref_seq_dict, ref_component_dict, ref_cluster_dict), ref_write_pickle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' import object from pickle\n')
            with open(input_file['ref_read_pickle'], 'rb') as ref_read_pickle:
                (ref_seq_dict, ref_component_dict, ref_cluster_dict) = pickle.load(ref_read_pickle)
    
        sys.stdout.write(Fore.MAGENTA + '[PROCESS]' + Style.RESET_ALL + ' generate report\n')
        ref_cluster_df_list = list()
        ref_seq_df = generate_report(input_file, ref_seq_dict, 'seq', 'ref')
        ref_component_df = generate_report(input_file, ref_component_dict, 'component', 'ref')
        for xprs_label in input_file['ref_xprs_label']:
            ref_cluster_df_list.append(generate_report(input_file, ref_cluster_dict, 'cluster', 'ref', xprs_label))
        ref_merge_df = merge_report(input_file, 'ref', ref_seq_df, ref_component_df, ref_cluster_df_list)
    
    if args.contig:
        sys.stdout.write(Fore.BLUE + '[CONTIG]\n')
        if 'contig_read_pickle' not in input_file:
            sys.stdout.write(Fore.BLUE + '[PROCESS]' + Style.RESET_ALL + ' object construction\n')
            contig_seq_dict = construct_sequence(input_file['contig_fasta'])
            #if 'contig_xprs_file' in input_file:
            sys.stdout.write(Fore.BLUE + '[PROCESS]' + Style.RESET_ALL + ' import expression value\n')
            read_expression(input_file, contig_seq_dict, 'contig')
            sys.stdout.write(Fore.BLUE + '[PROCESS]' + Style.RESET_ALL + ' calculation of sequence ambiguity\n')
            contig_self_blastn = filter_blastn(input_file['contig_blastn'])
            contig_self_match_dict = intersect_match(contig_self_blastn, contig_seq_dict, copy.deepcopy(contig_seq_dict))
            sys.stdout.write(Fore.BLUE + '[PROCESS]' + Style.RESET_ALL + ' construction of connected-component based graph\n')
            contig_uf, contig_component_dict = construct_graph(contig_seq_dict, contig_self_match_dict)
            sys.stdout.write(Fore.BLUE + '[PROCESS]' + Style.RESET_ALL + ' cluster expression\n')
            contig_cluster_dict = cluster_xprs(input_file, contig_seq_dict, contig_component_dict, 'contig')
            if 'contig_transrate' in input_file:
                sys.stdout.write(Fore.BLUE + '[PROCESS]' + Style.RESET_ALL + ' import transrate score\n')
                read_transrate(input_file, contig_seq_dict, 'contig')
            if 'contig_write_pickle' in input_file:   
                with open(input_file['contig_write_pickle'], 'wb') as contig_write_pickle:
                    pickle.dump((contig_seq_dict, contig_component_dict, contig_cluster_dict), contig_write_pickle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            sys.stdout.write(Fore.BLUE + '[PROCESS]' + Style.RESET_ALL + ' import object from pickle\n')
            with open(input_file['contig_read_pickle'], 'rb') as contig_read_pickle:
                (contig_seq_dict, contig_component_dict, contig_cluster_dict) = pickle.load(contig_read_pickle)
        
        sys.stdout.write(Fore.BLUE + '[PROCESS]' + Style.RESET_ALL + ' generate report\n')
        contig_cluster_df_list = list()
        contig_seq_df = generate_report(input_file, contig_seq_dict, 'seq', 'contig')
        contig_component_df = generate_report(input_file, contig_component_dict, 'component', 'contig')
        for xprs_label in input_file['contig_xprs_label']:
            contig_cluster_df_list.append(generate_report(input_file, contig_cluster_dict, 'cluster', 'contig', xprs_label))
        contig_merge_df = merge_report(input_file, 'contig', contig_seq_df, contig_component_df, contig_cluster_df_list)
        
    if args.reference and args.contig and args.match:
        sys.stdout.write(Fore.CYAN + '[MATCH]\n')
        match_blastn = filter_blastn(input_file['match_blastn'])
        sys.stdout.write(Fore.CYAN + '[PROCESS]' + Style.RESET_ALL + ' construction of sequence connection\n')
        match_dict = intersect_match(match_blastn, contig_seq_dict, ref_seq_dict)
        
        match_df = generate_report(input_file, match_dict, 'match')
        sys.stdout.write(Fore.CYAN + '[PROCESS]' + Style.RESET_ALL + ' merge report\n')
        merge_report(input_file, 'match', ref_merge_df=ref_merge_df, contig_merge_df=contig_merge_df, match_df=match_df)

    sys.stdout.write(Fore.GREEN + '[COMPLETE]' + Style.RESET_ALL + ' time elapse: ')
    sys.stdout.write(str(datetime.timedelta(seconds=round(time.time()-start_time))) + '\n')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('--gene', action='store_true', help='perform gene-isoform analysis on reference')
    parser.add_argument('--reference', action='store_true', help='perform reference analysis')
    parser.add_argument('--contig', action='store_true', help='perform contig analysis')
    parser.add_argument('--match', action='store_true', help='perform match analysis')    
    #parser.add_argument('--write', action='store_true', help='write reference pickle')
    #parser.add_argument('--read', action='store_true', help='read reference pickle')
    
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', type=str, help='input files in json format', required=True)
    
    args = parser.parse_args()
        
    main(args)