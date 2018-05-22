
import json
import pandas as pd
import numpy as np
import operator
import os
import re
import sys

from component import Component
from collections import defaultdict
from colorama import Fore, Style
from match import Match
from sequence import Sequence 
from sklearn.cluster import DBSCAN
from union_find import UnionFind

def process_and_validate_argument(args):
    terminate = False

    if args.match and not args.reference:
        sys.stderr.write(Fore.RED + '[WARNING]')
        sys.stderr.write(Style.RESET_ALL + ' ignore --match while --reference is not set\n')
    
    if args.match and not args.contig:
        sys.stderr.write(Fore.RED + '[WARNING]')
        sys.stderr.write(Style.RESET_ALL + ' ignore --match while --contig is not set\n')
    
    if os.path.exists(args.input):
        json_file = os.path.abspath(args.input)
        with open(json_file) as input_json:
            try:
                input_file = json.load(input_json)
            except:
                sys.stderr.write(Fore.RED + '[ERROR]' + Style.RESET_ALL + ' failed to load JSON format\n')
                terminate = True
            else:
                required_file = ['output_dir']
                
                if args.reference:
                    if 'ref_read_pickle' in input_file:
                        if 'ref_write_pickle' in input_file:
                            sys.stderr.write(Fore.RED + '[WARNING]')
                            sys.stderr.write(Style.RESET_ALL + ' ignore ref_write_pickle while ref_read_pickle is set\n')
                    else:
                        required_file.extend(['ref_fasta', 'ref_blastn', 'ref_xprs_file'])
                   
                if args.contig:
                    if 'contig_read_pickle' in input_file:
                        if 'contig_write_pickle' in input_file:
                            sys.stderr.write(Fore.RED + '[WARNING]')
                            sys.stderr.write(Style.RESET_ALL + ' ignore contig_write_pickle while contig_read_pickle is set\n')
                    else:
                        required_file.extend(['contig_fasta', 'contig_blastn', 'contig_xprs_file'])
                
                if args.match:
                    required_file.append('match_blastn')
                
                for target in ['ref', 'contig']:
                    if target + '_xprs_file' in required_file:
                        if target + '_xprs_file' not in input_file:
                            continue
                        else:
                            xprs_length = len(input_file[target + '_xprs_file'])
                        required_file.extend([target + '_xprs_file',
                                              target + '_xprs_header',
                                              target + '_xprs_name_col',
                                              target + '_xprs_value_col'])
                        if target + '_xprs_label' not in input_file:
                            input_file[target + '_xprs_label'] = list()
                            for i in range(0, xprs_length):
                                input_file[target + '_xprs_label'].append('xprs_' + '{0:02d}'.format(i+1))
                            sys.stderr.write(Fore.RED + '[WARNING]' + Style.RESET_ALL)
                            sys.stderr.write(' automatically fill ' + target + '_xprs_label with xprs_number\n')
                        if target + '_xprs_header' not in input_file:
                            input_file[target + '_xprs_header'] = [True] * xprs_length
                            sys.stderr.write(Fore.RED + '[WARNING]' + Style.RESET_ALL)
                            sys.stderr.write(' automatically fill ' + target + '_xprs_header with true\n')
                    elif target + '_xprs_file' in input_file:
                        sys.stderr.write(Fore.RED + '[WARNING]')
                        sys.stderr.write(Style.RESET_ALL + ' ignore ' + target + '_xprs_* while ' + target + '_read_pickle is set\n')
                    
                for required in required_file:
                    if required not in input_file:
                        sys.stderr.write(Fore.RED + '[ERROR]' + Style.RESET_ALL + ' required')
                        sys.stderr.write(Fore.RED + ' ' + required + Style.RESET_ALL + ' in input file\n')
                        terminate = True
                        
                for var in input_file.keys():
                    file = input_file[var]
                    
                    if var in ['ref_xprs_file', 'contig_xprs_file']:
                        for i in range(len(file)):
                            if os.path.exists(input_file[var][i]):
                                input_file[var][i] = os.path.abspath(input_file[var][i])
                            else:
                                sys.stderr.write(Fore.RED + '[ERROR] ' + var + ':' + str(i+1) + Style.RESET_ALL + ' does not exist\n') 
                                terminate = True
                        continue
                    elif isinstance(file, list):
                        continue
                        
                    if var in ['ref_write_pickle', 'contig_write_pickle', 'output_dir']:
                        path = os.path.dirname(file)
                        if var == 'output_dir':
                            path = os.path.abspath(file)
                        if not os.path.exists(path):
                            sys.stderr.write(Fore.RED + '[WARNING] ' + Style.RESET_ALL + 'create directory for ')
                            sys.stderr.write(Fore.RED + var + Style.RESET_ALL + '\n')
                            os.makedirs(path, exist_ok=True)
                        input_file[var] = os.path.abspath(file)
                        continue
                    
                    if os.path.exists(file):
                        input_file[var] = os.path.abspath(file)
                    else:
                        sys.stderr.write(Fore.RED + '[ERROR] ' + file + Style.RESET_ALL + ' does not exist\n') 
                        terminate = True     
    else:
        sys.stderr.write(Fore.RED + '[ERROR]' + Style.RESET_ALL + ' input JSON does not exsist\n')
        terminate = True

    if terminate:
        sys.exit(1)
    else:
        return input_file, args

def construct_sequence(file):
    length = 0
    sequences = dict()
    with open(file, 'r') as fasta:
        for i, line in enumerate(fasta):
            sline = line.strip()
            regex = re.match('>(\S+)', sline)
            if regex:
                if i != 0:
                    sequences[name] = Sequence(name, length)
                name = regex.group(1)
                length = 0
            else:        
                length += len(sline)
                
        sequences[name] = Sequence(name, length)
    
    return sequences

def gene_isoform_analysis(ref_gtf, seq_dict):

    gene_dict = dict()
    transcript_parent = dict()
    
    ref_gtf = pd.read_table(ref_gtf, sep='\t', header=None, low_memory=False)
    ref_gtf.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'header']
    
    for i, data in ref_gtf.iterrows():
        regex = re.match('transcript_id "(\S+)"; gene_id "(\S+)"', data['header'])
        if regex:
            transcript_name = regex.group(1)
            gene_name = regex.group(2)
        
            seq = seq_dict[transcript_name]
            seq.label['gene'] = gene_name
            
            gene = Component(gene_name)
            gene.add_member(seq)
            
            if gene_name in gene_dict:
                if seq not in gene_dict[gene_name].member:
                    gene_dict[gene_name].add_member(seq)
            else:
                gene_dict[gene_name] = gene
        
    return gene_dict

def filter_blastn(file, identity_t=70, evalue_t=1e-5, length_t=0):
    blastn = pd.read_table(file, sep='\t', header=None)
    blastn.columns = ['q_name', 'r_name', 'identity', 'm_length', 'mismatch', 'gap', 
                      'q_start', 'q_end', 'r_start', 'r_end', 'evalue', 'bitscore']
    length_f = blastn.loc[:, 'm_length'] >= length_t
    identity_f = blastn.loc[:, 'identity'] >= identity_t
    evalue_f = blastn.loc[:, 'evalue'] <= evalue_t
    duplicate_f = blastn.loc[:, 'q_name'] != blastn.loc[:, 'r_name']
    filtered_blastn = blastn.loc[length_f & identity_f & evalue_f & duplicate_f, :]
    
    return filtered_blastn

def intersect_match(blastn_dataframe, q_seq_dict, r_seq_dict):

    match_dict = dict()
    for i, data in blastn_dataframe.iterrows():
        
        q_seq = q_seq_dict[data['q_name']]
        r_seq = r_seq_dict[data['r_name']]
                
        match = Match(data, q_seq, r_seq)
        
        if match.m_name not in match_dict:
            match_dict[match.m_name] = match
        else:
            match_dict[match.m_name].extend(data, q_seq, r_seq)
        
    for m_name in match_dict.keys():
        match = match_dict[m_name]
        q_seq = q_seq_dict[match.q_name]
        r_seq = r_seq_dict[match.r_name]
        match.calculate_identity(q_seq, r_seq)
                   
    return match_dict

def construct_graph(seq_dict, match_dict, threshold=90):
    
    uf = UnionFind(seq_dict)
    component_dict = dict()
    
    for match in match_dict.values():
        q_seq = seq_dict[match.q_name]
        r_seq = seq_dict[match.r_name]
        if match.q_global_identity > threshold or match.r_global_identity > threshold:
            uf.union(q_seq.name, r_seq.name)
    
    uf.rename_component()
    
    for seq_name in seq_dict.keys():
        seq = seq_dict[seq_name]
        
        component_label = uf.component_label[seq_name]
        component_size = uf.component_size[component_label]
        
        seq.label['component'] = component_label
        
        component = Component(component_label)
        component.add_member(seq)
        
        if component_label in component_dict:
            component_dict[component_label].add_member(seq)
        else:
            component_dict[component_label] = component
          
    return uf, component_dict

def read_expression(input_file, seq_dict, target):

    for i, xprs_file in enumerate(input_file[target + '_xprs_file']):
        label = input_file[target + '_xprs_label'][i]
        name_col = int(input_file[target + '_xprs_name_col'][i]) - 1
        value_col = int(input_file[target + '_xprs_value_col'][i]) - 1
        header = input_file[target + '_xprs_header'][i]
    
        with open(xprs_file, 'r') as xprs_input:
            for j, data in enumerate(xprs_input):
                if header and j == 0:
                    continue
                col_data = data.split('\t')
                seq = seq_dict[col_data[name_col]]
                seq.xprs[label] = round(float(col_data[value_col]), 3)
        
    return
    
def read_transrate(input_file, seq_dict, target):
    transrate = pd.read_csv(input_file[target + '_transrate'])
    for i, data in transrate.iterrows():
        seq = seq_dict[data['contig_name']]
        seq.tr['good'] = round(data['p_good'], 3)
        seq.tr['bases_covered'] = round(data['p_bases_covered'], 3)
        seq.tr['seq_true'] = round(data['p_seq_true'], 3)
        seq.tr['score'] = round(data['score'], 3)
        seq.tr['not_segmented'] = round(data['p_not_segmented'], 3)
    return

def generate_xprs_summary(seq_dict, component_dict, gene_dict=None):
    for seq_name in seq_dict.keys():
        seq = seq_dict[seq_name]
        
        for label in seq.label:
            if label == 'component':
                component = component_dict[seq.label[label]]
            else:
                component = gene_dict[seq.label[label]]
            for xprs_label in seq.xprs:
                if component.get_total_xprs(xprs_label) != 0:
                    seq.contribute_xprs[xprs_label][label] = round(seq.xprs[xprs_label] / component.get_total_xprs(xprs_label), 3)
                    seq.relative_xprs[xprs_label][label] = round(seq.xprs[xprs_label] / component.get_maximum_xprs(xprs_label), 3)
                else:
                    seq.contribute_xprs[xprs_label][label] = 0 
                    seq.relative_xprs[xprs_label][label] = 0
    return

def generate_report(input_file, data_dict, data, target=None, xprs_label=None):
    target_list = ()
    table = list()
    header = list()
    if data == 'seq':
        for i, name in enumerate(sorted(data_dict.keys())):
            output_field = list()
            seq = data_dict[name]
            output_field.extend([seq.name, seq.length])
            if i == 0:
                header.extend([target + '_name', target + '_length'])
            
            for metric in sorted(seq.tr.keys()):
                output_field.append(seq.tr[metric])
                if i == 0:
                    header.append(target + '_tr_' + metric)
            
            for xprs_label in sorted(seq.xprs.keys()):
                output_field.append(seq.xprs[xprs_label])
                if i == 0:
                    header.append(target + '_xprs_' + xprs_label)
               
            for label in sorted(seq.label.keys(), reverse=True):
                output_field.append(seq.label[label])
                if i == 0:
                    header.append(target + '_' + label)
                for xprs_label in sorted(seq.xprs):
                    output_field.extend([seq.contribute_xprs[xprs_label][label], seq.relative_xprs[xprs_label][label]])
                    if i == 0:
                        header.extend([target + '_' + label + '_contribute_xprs_' + xprs_label,
                                       target + '_' + label + '_relative_xprs_' + xprs_label])
            
            table.append(output_field)
    elif data in ['component', 'gene']:
        for i, name in enumerate(sorted(data_dict.keys())):
            output_field = list()
            component = data_dict[name] 
            output_field.extend([component.name, len(component.member)])
            if i == 0:
                header.extend([target + '_' + data, target + '_' + data + '_size'])
            for label in sorted(component.member_xprs.keys()):
                output_field.extend([component.get_total_xprs(label),
                                     component.get_maximum_xprs(label),
                                     component.get_average_xprs(label)])
                if i == 0:
                    header.extend([target + '_' + data + '_tot_xprs_' + label,
                                   target + '_' + data + '_max_xprs_' + label,
                                   target + '_' + data + '_avg_xprs_' + label])
            table.append(output_field)
            
    elif data == 'match':
        for i, name in enumerate(sorted(data_dict.keys())):
            output_field = list()
            match = data_dict[name]
            output_field.extend([match.m_name, match.q_name, match.r_name,
                                 match.q_global_identity, match.r_global_identity])
            if i == 0:
                header.extend(['match_name', 'contig_name', 'ref_name', 'accuracy', 'recovery'])
            table.append(output_field)

    return pd.DataFrame(table, columns=header)

def merge_report(input_file, target,
                 seq_df=None, component_df=None, gene_df=None, 
                 ref_merge_df=None, contig_merge_df=None, match_df=None):
    
    if target in ['ref', 'contig']:
        if isinstance(gene_df, pd.DataFrame) and len(gene_df.index) != 0:
            merge_df = pd.merge(seq_df, gene_df, on=target+'_gene', how='left')
            merge_df = pd.merge(merge_df, component_df, on=target+'_component', how='left')
        else:
            merge_df = pd.merge(seq_df, component_df, on=target+'_component', how='left')
        merge_df.to_csv(input_file['output_dir'] + '/' + target + '.tsv', sep='\t', index=False)
    elif target == 'match':
        merge_df = pd.merge(match_df, contig_merge_df, on='contig_name', how='inner')
        merge_df = pd.merge(merge_df, ref_merge_df, on='ref_name', how='inner')
        merge_df.to_csv(input_file['output_dir'] + '/match.tsv', sep='\t', index=False)
    
    return merge_df