
import os
import pickle
import re

import numpy as np
import pandas as pd

from collections import defaultdict
from scipy.spatial.distance import cosine

class Sequence:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        
        self.expression_tpm = dict()
        
        self.ss = set()
        self.which_select = dict()
        self.select_which = dict()
        self.select_which_ordered = list()
        
        self.ss_depth = np.zeros(length).astype(float)
        
        self.select_which_depth = np.zeros(length).astype(float)
        self.which_select_depth = np.zeros(length).astype(float)

    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, target):
        return isinstance(target, Sequence) and self.name == target.name
    
class Match:
    def __init__(self, data, q_seq, r_seq):
        self.q_name = data['q_name']
        self.r_name = data['r_name']
        self.m_name = '.'.join([data['q_name'], data['r_name']])
        self.bitscore = data['bitscore']
        
        q_start = data['q_start']
        q_end = data['q_end']
        r_start = data['r_start']
        r_end = data['r_end']
        
        identity = data['identity'] / 100
        
        alignment = np.zeros(q_seq.length).astype(float)
        alignment[q_start-1:q_end] = identity
        self.q_m_depth = alignment
        
        alignment = np.zeros(r_seq.length).astype(float)
        alignment[r_start-1:r_end] = identity
        self.r_m_depth = alignment
        
    def __hash__(self):
        return hash(self.m_name)
        
    def __eq__(self, target):
        return isinstance(target, Match) and self.m_name == target.m_name
    
    def __gt__(self, target):
        return isinstance(target, Match) and np.sum(self.q_m_depth) > np.sum(target.q_m_depth)
    
    def extend(self, data, q_seq, r_seq):
        self.bitscore += data['bitscore']
        
        q_start = data['q_start']
        q_end = data['q_end']
        r_start = data['r_start']
        r_end = data['r_end']
        identity = data['identity'] / 100
        
        alignment = np.zeros(q_seq.length).astype(int)
        alignment[q_start-1:q_end] = identity
        self.q_m_depth = np.maximum(self.q_m_depth, alignment)
        
        alignment = np.zeros(r_seq.length).astype(int)
        alignment[r_start-1:r_end] = identity
        self.r_m_depth = np.maximum(self.r_m_depth, alignment)
        
        return

def find_match(blastn_dataframe, q_sequences, r_sequences):

    for i, data in blastn_dataframe.iterrows():
        
        q_name = data['q_name']
        r_name = data['r_name']
        q_seq = q_sequences[q_name]
        r_seq = r_sequences[r_name]
        
        q_start = data['q_start']
        r_start = data['r_start']
        q_end = data['q_end']
        r_end = data['r_end']
        
        identity = data['identity'] / 100
        bitscore = data['bitscore']
        
        alignment = np.zeros(q_seq.length).astype(int)
        alignment[q_start-1:q_end] = identity
        q_seq.select_which_depth = np.add(q_seq.select_which_depth, alignment)
        
        alignment = np.zeros(r_seq.length).astype(int)
        alignment[r_start-1:r_end] = identity
        r_seq.which_select_depth = np.add(r_seq.which_select_depth, alignment)
        
        match = Match(data, q_seq, r_seq)
        if q_name not in r_seq.which_select:
            r_seq.which_select[q_name] = match
        else:
            r_seq.which_select[q_name].extend(data, q_seq, r_seq)
        if r_name not in q_seq.select_which:
            q_seq.select_which[r_name] = match
        else:
            q_seq.select_which[r_name].extend(data, q_seq, r_seq)
        
    for q_name in q_sequences.keys():
        q_seq = q_sequences[q_name]
        if len(q_seq.select_which) != 0:
            ordered_matches = sorted(q_seq.select_which.values(), reverse=True)
            q_seq.select_which_ordered = ordered_matches
            
    return

def read_blastn(file, identity_t=95, evalue_t=1e-5, length_t=0):
    blastn = pd.read_table(file, sep='\t', header=None)
    blastn.columns = ['q_name', 'r_name', 'identity', 'm_length', 'mismatch', 'gap', 
                      'q_start', 'q_end', 'r_start', 'r_end', 'evalue', 'bitscore']
    length_f = blastn.loc[:, 'm_length'] >= length_t
    identity_f = blastn.loc[:, 'identity'] >= identity_t
    evalue_f = blastn.loc[:, 'evalue'] <= evalue_t
    duplicate_f = blastn.loc[:, 'q_name'] != blastn.loc[:, 'r_name']
    forward_f = blastn.loc[:, 'r_start'] <= blastn.loc[:, 'r_end']
    filtered_blastn = blastn.loc[length_f & identity_f & evalue_f & duplicate_f & forward_f, :]
    
    return filtered_blastn

def find_similar_subsequences(blastn_dataframe, sequences):    
    for i, data in blastn_dataframe.iterrows():
        
        q_seq = sequences[data['q_name']]
        r_seq = sequences[data['r_name']]
        q_start = data['q_start']
        q_end = data['q_end']
        identity = data['identity'] / 100
        
        alignment = np.zeros(q_seq.length).astype(int)
        alignment[q_start-1:q_end] = identity
        
        q_seq.ss_depth = np.add(q_seq.ss_depth, alignment)
        
        if r_seq.name not in q_seq.ss:
            q_seq.ss.add(r_seq.name)
        
    return

def construct_sequences(file):
    
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

def read_expression(sequences, assembly, base_dir):
    kallisto = pd.read_table(base_dir + '/' + assembly + "/kallisto/abundance.tsv", sep='\t')
    kallisto = kallisto.rename(columns={'target_id': 'name', 'tpm': 'kallisto'})
    kallisto = kallisto.loc[:, ('name', 'kallisto')]
    
    rsem = pd.read_table(base_dir + '/' + assembly + "/rsem/rsem.isoforms.results", sep='\t')
    rsem = rsem.rename(columns={'transcript_id': 'name', 'TPM': 'rsem'})
    rsem = rsem.loc[:, ('name', 'rsem')]
    
    salmon = pd.read_table(base_dir + '/' + assembly + "/salmon/quant.sf", sep='\t')
    salmon = salmon.rename(columns={'Name': 'name', 'TPM': 'salmon'})
    salmon = salmon.loc[:, ('name', 'salmon')]
                        
    expression_table = pd.merge(salmon, rsem, on='name', how='inner')
    expression_table = pd.merge(kallisto, expression_table, on='name', how='inner')
    
    if assembly == 'mRNA':
        flux_columns = ['locus', 'name', 'conding', 'length', 'molecular_fraction', 'molecular_count', 
                        'fragment_fraction', 'fragment_count', 'read_fraction', 'read_count', 'covered', 
                        'chi_square', 'variation']
        flux_simulator = pd.read_table(base_dir + "/simulation/flux_simulator.pro", sep='\t', header=None)
        flux_simulator.columns = flux_columns
        flux_simulator['read_per_nucleotide'] = flux_simulator['read_count'] / flux_simulator['length']
        flux_simulator['answer'] = 10 ** 6 * flux_simulator['read_per_nucleotide'] / np.sum(flux_simulator['read_per_nucleotide'])
        flux_simulator = flux_simulator.loc[:, ('name', 'answer')]
        
        expression_table = pd.merge(expression_table, flux_simulator, on='name', how='inner')

    for i, data in expression_table.iterrows():
        seq = sequences[data['name']]
        if assembly == 'mRNA':
            seq.expression_tpm['answer'] = data['answer']
        seq.expression_tpm['kallisto'] = data['kallisto']
        seq.expression_tpm['rsem'] = data['rsem']
        seq.expression_tpm['salmon'] = data['salmon']
        
    return

def generate_pickle(contigs, mRNAs, feature_dir):
        
    with open(feature_dir + '/contigs.pickle', 'wb') as output:
        pickle.dump(contigs, output, pickle.HIGHEST_PROTOCOL)
        
    with open(feature_dir + '/mRNAs.pickle', 'wb') as output:
        pickle.dump(mRNAs, output, pickle.HIGHEST_PROTOCOL)
    
    return

def generate_table(contigs, mRNAs, feature_dir):
    mRNA_list = list()
    t_match_list = list()
    for mRNA in mRNAs.values():
        if len(mRNA.select_which_ordered) != 0:
            for rank, match in enumerate(mRNA.select_which_ordered):
                t_match_list.append([
                        mRNA.name + ' aligned to ' + match.r_name,
                        mRNA.name,
                        match.r_name,
                        rank + 1,
                        np.sum(match.q_m_depth),
                        np.sum(match.r_m_depth),
                        match.q_m_depth[match.q_m_depth == 0].size,
                        match.r_m_depth[match.r_m_depth == 0].size])
                
        mRNA_list.append([
                mRNA.name,
                mRNA.length,
                mRNA.expression_tpm['answer'],
                len(mRNA.ss),
                np.sum(mRNA.ss_depth),
                mRNA.ss_depth[mRNA.ss_depth == 0].size,
                len(mRNA.select_which),
                np.sum(mRNA.select_which_depth),
                mRNA.select_which_depth[mRNA.select_which_depth == 0].size,
                mRNA.select_which_depth[mRNA.select_which_depth == 1].size,
                len(mRNA.which_select),
                np.sum(mRNA.which_select_depth),
                mRNA.which_select_depth[mRNA.which_select_depth == 0].size,
                mRNA.which_select_depth[mRNA.which_select_depth == 1].size])
    
    t_match_table = pd.DataFrame(t_match_list, columns=['m_name',
                                                        't_name',
                                                        'c_name',
                                                        't_rank',
                                                        'm_depth_sum_t',
                                                        'm_depth_sum_c',
                                                        'm_depth_zero_t',
                                                        'm_depth_zero_c'])
    
    mRNA_table = pd.DataFrame(mRNA_list, columns=['t_name', 
                                                  't_length',
                                                  't_answer_tpm',
                                                  't_ss_count',
                                                  't_ss_depth_sum',  
                                                  't_ss_depth_zero',
                                                  't_sw_count',
                                                  't_sw_depth_sum',
                                                  't_sw_depth_zero',
                                                  't_sw_depth_one',
                                                  't_ws_count',
                                                  't_ws_depth_sum',
                                                  't_ws_depth_zero',
                                                  't_ws_depth_one'])
    contig_list = list()
    c_match_list = list()
    for contig in contigs.values():
        if len(contig.select_which_ordered) != 0:
            for rank, match in enumerate(contig.select_which_ordered):
                c_match_list.append([
                        contig.name + ' aligned to ' + match.r_name,
                        match.r_name,
                        contig.name,
                        rank + 1,
                        np.sum(match.r_m_depth),
                        np.sum(match.q_m_depth),
                        match.r_m_depth[match.r_m_depth == 0].size,
                        match.q_m_depth[match.q_m_depth == 0].size])
                
        contig_list.append([
                contig.name,
                contig.length,
                contig.expression_tpm['kallisto'],
                contig.expression_tpm['rsem'],
                contig.expression_tpm['salmon'],
                len(contig.ss),
                np.sum(contig.ss_depth),
                contig.ss_depth[contig.ss_depth == 0].size,
                len(contig.select_which),
                np.sum(contig.select_which_depth),
                contig.select_which_depth[contig.select_which_depth == 0].size,
                contig.select_which_depth[contig.select_which_depth == 1].size,
                len(contig.which_select),
                np.sum(contig.which_select_depth),
                contig.which_select_depth[contig.which_select_depth == 0].size,
                contig.which_select_depth[contig.which_select_depth == 1].size])
    
    c_match_table = pd.DataFrame(c_match_list, columns=['m_name',
                                                        't_name',
                                                        'c_name',
                                                        'c_rank',
                                                        'm_depth_sum_t',
                                                        'm_depth_sum_c',
                                                        'm_depth_zero_t',
                                                        'm_depth_zero_c'])
    
    contig_table = pd.DataFrame(contig_list, columns=['c_name', 
                                                      'c_length',
                                                      'c_kallisto_tpm',
                                                      'c_rsem_tpm',
                                                      'c_salmon_tpm',
                                                      'c_ss_count',
                                                      'c_ss_depth_sum',
                                                      'c_ss_depth_zero',
                                                      'c_sw_count',
                                                      'c_sw_depth_sum',
                                                      'c_sw_depth_zero',
                                                      'c_sw_depth_one',
                                                      'c_ws_count',
                                                      'c_ws_depth_sum',
                                                      'c_ws_depth_zero',
                                                      'c_ws_depth_one'])
    
    t_feature_table = pd.merge(t_match_table, mRNA_table, on='t_name', how='left')
    t_feature_table = pd.merge(t_feature_table, contig_table, on='c_name', how='left')
    
    c_feature_table = pd.merge(c_match_table, mRNA_table, on='t_name', how='left')
    c_feature_table = pd.merge(c_feature_table, contig_table, on='c_name', how='left')
    
    mRNA_table.to_csv(feature_dir + '/mRNA.tsv', sep='\t', index=False)
    contig_table.to_csv(feature_dir + '/contig.tsv', sep='\t', index=False)
    
    t_match_table.to_csv(feature_dir + '/t_match.tsv', sep='\t', index=False)
    c_match_table.to_csv(feature_dir + '/c_match.tsv', sep='\t', index=False)
    
    t_feature_table.to_csv(feature_dir + '/t_feature.tsv', sep='\t', index=False)
    c_feature_table.to_csv(feature_dir + '/c_feature.tsv', sep='\t', index=False)
    
    return t_feature_table, c_feature_table

if True:
    
    species = 'yeast_50x'
    assembly = 'trinity'
    
    base_dir = "/home/dn070017/projects/QuantEval/simulation/" + species
    mRNA_dir = base_dir + "/mRNA/"
    contig_dir = base_dir + "/" + assembly + "/"
    feature_dir = contig_dir + '/features/'
    
    if not os.path.exists(feature_dir):
        os.makedirs(feature_dir)
    
    contigs = contig_dir + "trinity.fasta"
    contigs = construct_sequences(contigs)
    
    mRNAs = mRNA_dir + "/mRNA.fasta"
    mRNAs = construct_sequences(mRNAs)
    
    read_expression(mRNAs, 'mRNA', base_dir)
    read_expression(contigs, assembly, base_dir)
    
    c_self_blastn = contig_dir + "/blastn/self.tsv"
    c_self_blastn = read_blastn(c_self_blastn)
    
    m_self_blastn = mRNA_dir + "/blastn/self.tsv"
    m_self_blastn = read_blastn(m_self_blastn)
    
    m_c_blastn = contig_dir + "/blastn/mRNA_to_contig.tsv"
    m_c_blastn = read_blastn(m_c_blastn)
    
    c_m_blastn = contig_dir + "/blastn/contig_to_mRNA.tsv"
    c_m_blastn = read_blastn(c_m_blastn)
    
    find_similar_subsequences(c_self_blastn, contigs)
    find_similar_subsequences(m_self_blastn, mRNAs)
    
    find_match(m_c_blastn, mRNAs, contigs)
    find_match(c_m_blastn, contigs, mRNAs)
    
    generate_pickle(contigs, mRNAs, feature_dir)
    
    t_feature_table, c_feature_table = generate_table(contigs, mRNAs, feature_dir)
    
    print('done')