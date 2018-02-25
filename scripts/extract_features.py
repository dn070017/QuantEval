#!/usr/bin/env python3

import copy
import os
import pickle
import re
import sys
import time

import numpy as np
import pandas as pd

from collections import defaultdict
from scipy.spatial.distance import cosine

class UnionFind:
    
    def __init__(self):
        self.parent = dict()
        self.size = dict()
    
    def find(self, a):
        if a not in self.parent:
            self.parent[a] = a
            self.size[a] = 1
        while a != self.parent[a]:
            a = self.parent[a]
        return a
    
    def union(self, a, b):
        a_root = self.find(a)
        b_root = self.find(b)
        if a_root == b_root:
            return
        else:
            if self.size[b_root] > self.size[a_root]:
                self.parent[a_root] = b_root
                self.size[b_root] += self.size[a_root]
            else:
                self.parent[b_root] = a_root
                self.size[a_root] += self.size[b_root]
         
    def get_size(self, a):
        return self.size[self.find(a)]
    
    def root(self):
        i = 0
        roots = dict()
        for k, v in self.parent.items():
            if k == v:
                roots[k] = i
                i += 1
        return roots

class Sequence:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        
        self.cc_label = -1
        self.cc_size = 1
        
        self.gene_name = np.nan
        self.gene_cluster_size = np.nan
        
        self.expression_tpm = dict()
        
        self.match_seq = dict()
        self.match_seq['ss'] = set()
        self.match_seq['tc'] = set()
        
        self.match_ordered = dict()
        self.match_ordered['ss'] = list()
        self.match_ordered['tc'] = list()
        
        self.match_count = dict()
        self.match_count['ss'] = np.zeros(length).astype(int)
        self.match_count['tc'] = np.zeros(length).astype(int)
        
        self.match_depth = dict()
        self.match_depth['ss'] = np.zeros(length).astype(float)
        self.match_depth['tc'] = np.zeros(length).astype(float)
        
        self.tr_good = np.nan
        self.tr_bases_covered = np.nan 
        self.tr_seq_true = np.nan 
        self.tr_score =np.nan 
        self.tr_not_segmented = np.nan
        
        self.corset_label = np.nan
        self.corset_size = np.nan

    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, target):
        return isinstance(target, Sequence) and self.name == target.name
    
class Match:
    def __init__(self, data, q_seq, r_seq):
        self.q_name = data['q_name']
        self.r_name = data['r_name']
        self.m_name = '.'.join([data['q_name'], data['r_name']])
        self.q_rank = np.nan
        self.r_rank = np.nan
        
        q_start = data['q_start']
        q_end = data['q_end']
        r_start = data['r_start']
        r_end = data['r_end']
        
        if r_start > r_end:
            (r_start, r_end) = (r_end, r_start)
        
        identity = data['identity'] / 100
        
        alignment = np.zeros(q_seq.length).astype(float)
        alignment[q_start-1:q_end] = identity
        self.q_depth = alignment
        
        count = np.zeros(q_seq.length).astype(int)
        count[q_start-1:q_end] = 1
        self.q_count = count
        
        alignment = np.zeros(r_seq.length).astype(float)
        alignment[r_start-1:r_end] = identity
        self.r_depth = alignment
        
        count = np.zeros(r_seq.length).astype(int)
        count[r_start-1:r_end] = 1
        self.r_count = count
        
    def __hash__(self):
        return hash(self.m_name)
        
    def __eq__(self, target):
        return isinstance(target, Match) and self.m_name == target.m_name
    
    def __gt__(self, target):
        return isinstance(target, Match) and np.sum(self.q_depth) > np.sum(target.q_depth)
    
    def extend(self, data, q_seq, r_seq):
        
        q_start = data['q_start']
        q_end = data['q_end']
        r_start = data['r_start']
        r_end = data['r_end']
        identity = data['identity'] / 100
        
        if r_start > r_end:
            (r_start, r_end) = (r_end, r_start)
        
        alignment = np.zeros(q_seq.length).astype(float)
        alignment[q_start-1:q_end] = identity
        self.q_depth = np.maximum(self.q_depth, alignment)
        
        count = np.zeros(q_seq.length).astype(int)
        count[q_start-1:q_end] = 1
        self.q_count = np.maximum(self.q_count, count)
        
        alignment = np.zeros(r_seq.length).astype(float)
        alignment[r_start-1:r_end] = identity
        self.r_depth = np.maximum(self.r_depth, alignment)
        
        count = np.zeros(r_seq.length).astype(int)
        count[r_start-1:r_end] = 1
        self.r_count = np.maximum(self.r_count, count)
        
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

def find_match(blastn_dataframe, q_sequences, r_sequences, mode):

    matches = dict()
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
        
        match = Match(data, q_seq, r_seq)
        
        if match.m_name not in matches:
            matches[match.m_name] = match
        else:
            matches[match.m_name].extend(data, q_seq, r_seq)
        
    for m_name in matches.keys():
        match = matches[m_name]
        q_name = match.q_name
        r_name = match.r_name
        q_seq = q_sequences[q_name]
        r_seq = r_sequences[r_name]
        
        q_seq.match_seq[mode].add(r_name)
        q_seq.match_ordered[mode].append(match)
        q_seq.match_depth[mode] = np.add(q_seq.match_depth[mode], match.q_depth)
        q_seq.match_count[mode] = np.add(q_seq.match_count[mode], match.q_count)
            
        r_seq.match_seq[mode].add(q_name)
        r_seq.match_ordered[mode].append(match)
        r_seq.match_depth[mode] = np.add(r_seq.match_depth[mode], match.r_depth)
        r_seq.match_count[mode] = np.add(r_seq.match_count[mode], match.r_count)
          
    for q_name in q_sequences.keys():
        q_seq = q_sequences[q_name]
        if len(q_seq.match_seq[mode]) != 0:
            q_seq.match_ordered[mode].sort(reverse=True)
            for i, match in enumerate(q_seq.match_ordered[mode]):
                if match.q_name == q_name:
                    match.q_rank = i + 1
                    
    for r_name in r_sequences.keys():
        r_seq = r_sequences[r_name]
        if len(r_seq.match_seq[mode]) != 0:
            r_seq.match_ordered[mode].sort(reverse=True)
            for i, match in enumerate(r_seq.match_ordered[mode]):
                match.r_rank = i + 1
           
    return matches

def read_blastn(file, identity_t=70, evalue_t=1e-5, length_t=0):
    blastn = pd.read_table(file, sep='\t', header=None)
    blastn.columns = ['q_name', 'r_name', 'identity', 'm_length', 'mismatch', 'gap', 
                      'q_start', 'q_end', 'r_start', 'r_end', 'evalue', 'bitscore']
    length_f = blastn.loc[:, 'm_length'] >= length_t
    identity_f = blastn.loc[:, 'identity'] >= identity_t
    evalue_f = blastn.loc[:, 'evalue'] <= evalue_t
    duplicate_f = blastn.loc[:, 'q_name'] != blastn.loc[:, 'r_name']
    filtered_blastn = blastn.loc[length_f & identity_f & evalue_f & duplicate_f, :]
    
    return filtered_blastn

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
        else:
            seq.expression_tpm['answer'] = np.nan
        seq.expression_tpm['kallisto'] = data['kallisto']
        seq.expression_tpm['rsem'] = data['rsem']
        seq.expression_tpm['salmon'] = data['salmon']
        
    return

def read_transrate(sequences, assembly, base_dir):
    transrate = pd.read_csv(base_dir + '/' + assembly + "/transrate/" + assembly + "/contigs.csv")
    for i, data in transrate.iterrows():
        sequences[data['contig_name']].tr_good = data['p_good']
        sequences[data['contig_name']].tr_bases_covered = data['p_bases_covered']
        sequences[data['contig_name']].tr_seq_true = data['p_seq_true']
        sequences[data['contig_name']].tr_score = data['score']
        sequences[data['contig_name']].tr_not_segmented = data['p_not_segmented']
    return

def read_corset(sequences, assembly, base_dir):
    corset_clusters = pd.read_table(base_dir + '/' + assembly + '/corset/corset-clusters.txt', sep='\t', header=None)
    corset_clusters.columns = ['name', 'cluster']
    count = corset_clusters['cluster'].value_counts()
    for i, data in corset_clusters.iterrows():
        cluster = data['cluster']
        cluster_size = count[cluster]
        sequences[data['name']].corset_label = cluster
        sequences[data['name']].corset_size = cluster_size
        
    i = 0
    for name in sequences.keys():
        seq = sequences[name]
        if np.isnan(seq.corset_size):
            seq.corset_label = 'NA-Cluster-' + str(i) + '.0'
            seq.corset_size = 1
            i += 1
    return

def match_gene_isoform(ref_dir, sequences):
    gtf = ref_dir + '/flux_simulator_clean.gtf'

    transcript_count = set()
    isoform_count = dict()
    transcript_parent = dict()
    gtf = pd.read_table(gtf, sep='\t', header=None, low_memory=False)
    gtf.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'header']
    for i, data in gtf.iterrows():
        regex = re.match('transcript_id "(\S+)"; gene_id "(\S+)"', data['header'])
        if regex:
            transcript_parent[regex.group(1)] = regex.group(2)
            if regex.group(1) in transcript_count:
                continue
            else:
                transcript_count.add(regex.group(1))
                if regex.group(2) in isoform_count:
                    isoform_count[regex.group(2)] += 1
                else:
                    isoform_count[regex.group(2)] = 1
                
    for name in sequences.keys():
        seq = sequences[name]
        seq.gene_name = transcript_parent[name]
        seq.gene_cluster_size = isoform_count[transcript_parent[name]]
        
    return

def return_match(match):
    return([match.q_name + ' aligned to ' + match.r_name, match.q_name, match.r_name,
            match.q_rank, match.r_rank, np.around(np.sum(match.q_depth), 2), np.around(np.sum(match.r_depth), 2),
            match.q_depth[match.q_depth == 0].size, match.r_depth[match.r_depth == 0].size])

def return_seq(seq, seq_type):
    data = [seq.name, seq.length, np.around(seq.expression_tpm['answer'], 2),
            np.around(seq.expression_tpm['kallisto'], 2), np.around(seq.expression_tpm['rsem'], 2),
            np.around(seq.expression_tpm['salmon'], 2), len(seq.match_seq['ss']), 
            np.around(np.sum(seq.match_depth['ss']), 2), seq.match_count['ss'][seq.match_count['ss'] == 0].size,
            len(seq.match_seq['tc']), np.around(np.sum(seq.match_depth['tc']), 2),
            seq.match_count['tc'][seq.match_count['tc'] == 0].size,
            seq.match_count['tc'][seq.match_count['tc'] == 1].size,
            seq.tr_good, seq.tr_bases_covered, seq.tr_seq_true, seq.tr_score, seq.tr_not_segmented,
            seq.corset_label, seq.corset_size]
    if seq_type == 'mRNA':
        return(data + [seq.gene_name, seq.gene_cluster_size])
    else:
        return(data)

def generate_table(contigs, mRNAs, tc_matches, feature_dir):
    
    t_ss_list = list()
    mRNA_list = list()
    for mRNA in mRNAs.values():
        if len(mRNA.match_seq['ss']) != 0:
            for match in mRNA.match_ordered['ss']:
                t_ss_list.append(return_match(match))
        mRNA_list.append(return_seq(mRNA, 'mRNA'))
    
    c_ss_list = list()
    contig_list = list()
    for contig in contigs.values():
        if len(contig.match_seq['ss']) != 0:
            for match in contig.match_ordered['ss']:
                c_ss_list.append(return_match(match))
        contig_list.append(return_seq(contig, 'contig'))
        
    tc_list = list()
    for match in tc_matches.values():
        tc_list.append(return_match(match))
        
    t_ss_table = pd.DataFrame(t_ss_list, columns=['m_name', 't_name', 't_name_r', 't_ss_rank', 't_ss_rank_r',
                                                  't_ss_m_depth_sum', 't_ss_m_depth_sum_r', 't_ss_m_depth_zero',
                                                  't_ss_m_depth_zero_r'])
    
    c_ss_table = pd.DataFrame(c_ss_list, columns=['m_name', 'c_name', 'c_name_r', 'c_ss_rank', 'c_ss_rank_r',
                                                  'c_ss_m_depth_sum', 'c_ss_m_depth_sum_r', 'c_ss_m_depth_zero',
                                                  'c_ss_m_depth_zero_r'])
    
    tc_table = pd.DataFrame(tc_list, columns=['m_name', 'c_name', 't_name', 'c_tc_rank', 't_tc_rank',
                                              'c_tc_m_depth_sum', 't_tc_m_depth_sum', 'c_tc_m_depth_zero',
                                              't_tc_m_depth_zero'])

    mRNA_table = pd.DataFrame(mRNA_list, columns=['t_name', 't_length', 't_answer_tpm', 't_kallisto_tpm',
                                                  't_rsem_tpm', 't_salmon_tpm', 't_ss_count', 't_ss_depth_sum',  
                                                  't_ss_depth_zero', 't_tc_count', 't_tc_depth_sum', 
                                                  't_tc_count_zero', 't_tc_count_one', 't_tr_good',
                                                  't_tr_bases_covered', 't_tr_seq_true', 't_tr_score',
                                                  't_tr_not_segmented', 't_corset_label', 't_corset_size',
                                                  't_gene_name', 't_gene_cluster_size'])
    
    contig_table = pd.DataFrame(contig_list, columns=['c_name', 'c_length', 'c_answer_tpm', 'c_kallisto_tpm',
                                                      'c_rsem_tpm', 'c_salmon_tpm', 'c_ss_count', 'c_ss_depth_sum',  
                                                      'c_ss_depth_zero', 'c_tc_count', 'c_tc_depth_sum',
                                                      'c_tc_count_zero', 'c_tc_count_one', 'c_tr_good', 
                                                      'c_tr_bases_covered', 'c_tr_seq_true', 'c_tr_score', 
                                                      'c_tr_not_segmented', 'c_corset_label', 'c_corset_size'])
    contig_table = contig_table.drop(['c_answer_tpm'], axis=1)
    
    feature_table = pd.merge(tc_table, mRNA_table, on='t_name', how='left')
    feature_table = pd.merge(feature_table, contig_table, on='c_name', how='left')
    
    mRNA_table.to_csv(feature_dir + '/mRNA.tsv', sep='\t', index=False)
    contig_table.to_csv(feature_dir + '/contig.tsv', sep='\t', index=False)
    
    t_ss_table.to_csv(feature_dir + '/mRNA_ss.tsv', sep='\t', index=False)
    c_ss_table.to_csv(feature_dir + '/contig_ss.tsv', sep='\t', index=False)
    tc_table.to_csv(feature_dir + '/CTPs.tsv', sep='\t', index=False)
    
    feature_table.to_csv(feature_dir + '/features.tsv', sep='\t', index=False)
    
    return feature_table, mRNA_table, contig_table, t_ss_table, c_ss_table, tc_table

def generate_mRNA_table(mRNAs, feature_dir):
    mRNA_list = list()
    for mRNA in mRNAs.values():                
        mRNA_list.append([
                mRNA.name,
                mRNA.length,
                np.around(mRNA.expression_tpm['answer'], 2),
                np.around(mRNA.expression_tpm['kallisto'], 2),
                np.around(mRNA.expression_tpm['rsem'], 2),
                np.around(mRNA.expression_tpm['salmon'], 2),
                len(mRNA.match_seq['ss']),
                np.around(np.sum(mRNA.match_depth['ss']), 2),
                mRNA.match_count['ss'][mRNA.match_count['ss'] == 0]])
    
    mRNA_table = pd.DataFrame(mRNA_list, columns=['t_name', 
                                                  't_length',
                                                  't_answer_tpm',
                                                  't_kallisto_tpm',
                                                  't_rsem_tpm',
                                                  't_salmon_tpm',
                                                  't_ss_count',
                                                  't_ss_depth_sum',  
                                                  't_ss_depth_zero'])
    
    mRNA_table.to_csv(feature_dir + '/mRNA.tsv', sep='\t', index=False)
   
    return mRNA_table

def main(base_dir, ref_dir, assembly):
#if True:
    print('start main program')
    start_time = time.time()

    #species = 'mouse_50x'
    #assembly = 'rnaspades'
    #assembly = 'transabyss'
    #assembly = 'trinity'
    
    #base_dir = "/home/dn070017/projects/QuantEval/simulation/dog_50x/"
    #ref_dir = "/home/dn070017/projects/QuantEval/reference/dog/"
    mRNA_dir = base_dir + "/mRNA/"
    contig_dir = base_dir + "/" + assembly + "/"
    
    mRNAs = mRNA_dir + "/mRNA.fasta"
    mRNAs = construct_sequences(mRNAs)
    read_expression(mRNAs, 'mRNA', base_dir)
    read_transrate(mRNAs, 'mRNA', base_dir)
    read_corset(mRNAs, 'mRNA', base_dir)
    
    print('    - start analyzing gene-isoform relation...')
    match_gene_isoform(ref_dir, mRNAs)
    
    print('    - start analyzing similar subsequences in transcripts...')
    m_self_blastn = mRNA_dir + "/blastn/self.tsv"
    m_self_blastn = read_blastn(m_self_blastn)
    t_ss_matches = find_match(m_self_blastn, mRNAs, copy.deepcopy(mRNAs), 'ss')
    
    if assembly == 'mRNA':
        feature_dir = mRNA_dir + '/features/'
        if not os.path.exists(feature_dir):
            os.makedirs(feature_dir)
        mRNA_table = generate_mRNA_table(mRNAs, feature_dir)
    
    else:
        feature_dir = contig_dir + '/features/'
        if not os.path.exists(feature_dir):
            os.makedirs(feature_dir)
        
        contigs = contig_dir + '/' + assembly + ".fasta"
        contigs = construct_sequences(contigs)
        read_expression(contigs, assembly, base_dir)
        read_transrate(contigs, assembly, base_dir)
        read_corset(contigs, assembly, base_dir)
        
        print('    - start analyzing similar subsequences in contigs...')
        c_self_blastn = contig_dir + "/blastn/self.tsv"
        c_self_blastn = read_blastn(c_self_blastn)
        c_ss_matches = find_match(c_self_blastn, contigs, copy.deepcopy(contigs), 'ss')
        
        print('    - start analyzing matches between transcripts and contigs...')
        c_m_blastn = contig_dir + "/blastn/contig_to_mRNA.tsv"
        c_m_blastn = read_blastn(c_m_blastn)
        tc_matches = find_match(c_m_blastn, contigs, mRNAs, 'tc')
        
        print('    - start generating summary tables...')        
        feature_table, mRNA_table, contig_table, t_ss_table, c_ss_table, tc_table = generate_table(contigs, mRNAs, tc_matches, feature_dir)
    
    print('finished, time elapsed: ', np.around(time.time() - start_time, 3), ' seconds', sep='')
    
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('usage: ', os.path.basename(__file__), ' [base directory] [reference directory] [assembly]', sep='')
        sys.exit(1)
    base_dir = sys.argv[1]
    ref_dir = sys.argv[2]
    assembly = sys.argv[3]
    main(base_dir, ref_dir, assembly)
    
# DEPRECATED FUNCTION #
def generate_pickle(contigs, mRNAs, feature_dir):
    with open(feature_dir + '/dev_contigs.pickle', 'wb') as output:
        pickle.dump(contigs, output, pickle.HIGHEST_PROTOCOL)
        
    with open(feature_dir + '/dev_mRNAs.pickle', 'wb') as output:
        pickle.dump(mRNAs, output, pickle.HIGHEST_PROTOCOL)
    
    return

# DEPRECATED FUNCTION #
def group_cc(seqs, matches):
    uf = UnionFind()
    for match in matches.values():
        if match.q_depth[match.q_depth>=0.95].size/match.q_count.size > 0.99 and \
           match.r_depth[match.r_depth>=0.95].size/match.r_count.size > 0.99:
            uf.union(match.q_name, match.r_name)
    
    root = uf.root()
    for cc in uf.parent.keys():
        root_name = uf.find(cc)
        seqs[cc].cc_label = root[root_name]
        seqs[cc].cc_size = uf.get_size(root_name)

    return