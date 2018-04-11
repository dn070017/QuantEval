#!/usr/bin/env python3

import copy
import operator
import os
import pickle
import re
import sys
import time

import numpy as np
import pandas as pd

from collections import defaultdict
from scipy.spatial.distance import cosine
from sklearn.cluster import DBSCAN

class UnionFind:
   
    def __init__(self, seqs):
        self.parent = dict()
        self.size = dict()
        for name in seqs.keys():
            self.parent[name] = name
            self.size[name] = 1
   
    def find(self, a):
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
    
    def get_root_label(self):
        i = 1
        root_label = dict()
        for k, v in self.parent.items():
            if k == v:
                root_label[k] = 'cluster_' + str(i)
                i += 1
        return root_label

class Sequence:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.xprs = dict()
        self.clustered = False
        self.cluster_label = 'unclustered'
        self.cluster_size = 0
        self.match_seq = set()
        self.match_ordered = list()
        self.match_count = np.zeros(length).astype(int)
        self.match_depth = np.zeros(length).astype(float)

    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, target):
        return isinstance(target, Sequence) and self.name == target.name
    
class Match:
    def __init__(self, data, q_seq, r_seq):
        self.q_name = data['q_name']
        self.r_name = data['r_name']
        
        q_start = data['q_start']
        q_end = data['q_end']
        r_start = data['r_start']
        r_end = data['r_end']
        
        if r_start > r_end:
            (r_start, r_end) = (r_end, r_start)
            self.orientation = 'R'
        else:
            self.orientation = 'F'
    
        self.m_name = '.'.join([self.q_name, self.r_name, self.orientation])
    
        identity = data['identity'] / 100
    
        alignment = np.zeros(q_seq.length).astype(float)
        alignment[q_start-1:q_end] = identity
        self.q_depth = alignment
        
        alignment = np.zeros(r_seq.length).astype(float)
        alignment[r_start-1:r_end] = identity
        self.r_depth = alignment

    def __hash__(self):
        return hash(self.m_name)
        
    def __eq__(self, target):
        return isinstance(target, Match) and self.m_name == target.m_name
    
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
        
        alignment = np.zeros(r_seq.length).astype(float)
        alignment[r_start-1:r_end] = identity
        self.r_depth = np.maximum(self.r_depth, alignment)

        return
    
    def calculate_identity(self, q_seq, r_seq):
        self.q_global_identity = np.sum(self.q_depth) / q_seq.length
        self.r_global_identity = np.sum(self.r_depth) / r_seq.length
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

def find_match(blastn_dataframe, q_sequences, r_sequences):

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
        match.calculate_identity(q_seq, r_seq)
           
    return matches

def read_blastn(file, identity_t, evalue_t=1e-5, length_t=0):
    blastn = pd.read_table(file, sep='\t', header=None)
    blastn.columns = ['q_name', 'r_name', 'identity', 'm_length', 'mismatch', 'gap', 
                      'q_start', 'q_end', 'r_start', 'r_end', 'evalue', 'bitscore']
    length_f = blastn.loc[:, 'm_length'] >= length_t
    identity_f = blastn.loc[:, 'identity'] >= identity_t * 100
    evalue_f = blastn.loc[:, 'evalue'] <= evalue_t
    duplicate_f = blastn.loc[:, 'q_name'] != blastn.loc[:, 'r_name']
    filtered_blastn = blastn.loc[length_f & identity_f & evalue_f & duplicate_f, :]
    
    return filtered_blastn

def read_xprs(contigs, xprs, n_col, x_col, header=True, sep='\t'):
    if header:
        header=0
    else:
        header=None
    xprs = pd.read_csv(xprs, header=header, sep=sep)
    xprs = xprs.set_index(xprs.iloc[:,n_col])
    xprs = xprs.to_dict()[xprs.columns[x_col]]
    for contig in contigs.values():
        contig.xprs = xprs[contig.name]
    return

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

def connected_component(seqs, matches, threshold):
    uf = UnionFind(seqs)
    clusters = defaultdict(list)
    
    for match in matches.values():
        q_seq = seqs[match.q_name]
        r_seq = seqs[match.r_name]
        if match.q_global_identity > threshold or match.r_global_identity > threshold:
            uf.union(q_seq.name, r_seq.name)
    
    root_label = uf.get_root_label()
    for seq_name in uf.parent.keys():
        root_name = uf.find(seq_name)
        seqs[seq_name].cluster_label = root_label[root_name]
        seqs[seq_name].cluster_size = uf.get_size(root_name)
        clusters[root_label[root_name]].append(seqs[seq_name])

    return uf, clusters

def xprs_distance(cluster, xprs):
    length = len(cluster)
    if np.max(np.abs(xprs) != 0):
        norm_xprs = xprs / np.sum(np.abs(xprs))
    else:
        norm_xprs = xprs
    dist = np.ones((length, length))
    for i in range(0, length):
        for j in range(0, length):
            if i == j:
                dist[i, j] = 0
            else:
                q = cluster[i]
                r = cluster[j]
                
                match_score = np.absolute(norm_xprs[i] - norm_xprs[j])
        
                dist[i, j] = match_score  
    return(dist)

def quant_eval(clusters, dbscan_eps, out_file):
    
    output_list = list()
    
    for cluster_name, cluster in clusters.items():

        cluster = np.array(cluster)
        xprs = list()
        for seq in cluster:
            xprs.append(seq.xprs)

        xprs = np.array(xprs)  
        dist = xprs_distance(cluster, xprs)

        if len(xprs) <= 2:
            min_samples = 2
        else:
            min_samples = 3
        
        dbscan = DBSCAN(metric="precomputed", eps=dbscan_eps, min_samples=min_samples).fit(dist)

        outlier_id = 1
        labels = list()
        for label in dbscan.labels_:
            if label == -1:
                labels.append(-1 * outlier_id)
                outlier_id += 1
            else:
                labels.append(label)

        total_sums = 0
        cluster_sums = dict()
        cluster_means = dict()
        cluster_counts = dict()
        cluster_rank = dict()
        for label in np.unique(labels):
            indices = np.where(labels==label)
            total_sums += np.sum(np.take(xprs, indices))
            cluster_sums[label] = np.sum(np.take(xprs, indices))
            cluster_means[label] = np.mean(np.take(xprs, indices))
            cluster_counts[label] = len(indices[0])

        sorted_means = sorted(cluster_means.items(), key=operator.itemgetter(1), reverse=True)    
        rank = 1
        for label, means in sorted_means:
            cluster_rank[label] = rank
            rank += 1 
        
        for i, seq in enumerate(cluster):
            label = labels[i]
            count = cluster_counts[label]
            sums = cluster_sums[label]
            means = cluster_means[label]
            rank = cluster_rank[label]
            tmp_list = [seq.name, cluster_name, seq.cluster_size, label, rank, count, seq.xprs, total_sums, sums]

            if rank == 1:
                if count == 1:
                    tmp_list.extend([seq.xprs, 'P'])
                else:
                    tmp_list.extend([sums, 'U'])
            else:
                tmp_list.extend([seq.xprs, 'S'])
            
            output_list.append(tmp_list)
        
    header = ['c_name', 'sequence_cluster', 'size_of_sequence_cluster', 'xprs_cluster', 'rank_of_xprs_cluster', 
              'size_of_xprs_cluster', 'original_xprs', 'sequence_cluster_xprs', 'xprs_cluster_xprs',  
              'corrected_xprs', 'category']

    output_pd = pd.DataFrame(output_list, columns=header)
    output_pd.to_csv(out_file, sep='\t', index=False)
    return

def main(blastn, fasta, xprs, n_col, x_col, header, sep, blast_identity, cluster_identity, dbscan_eps, out_file):
#if True:
    print('start main program')
    start_time = time.time()

    #species = 'dog_50x'
    #assembly = 'rnaspades'
    #assembly = 'transabyss'
    #assembly = 'trinity'
    
    #base_dir = "/home/dn070017/projects/QuantEval/simulation/dog_50x/"
    #contig_dir = base_dir + "/" + assembly + "/"
    #out_dir = contig_dir + "/quateval"
    #xprs = contig_dir + '/kallisto/abundance.tsv'
    #fasta = contig_dir + '/' + assembly + ".fasta"
    #blastn = contig_dir + "/blastn/self.tsv"
    
    #n_col = 0
    #x_col = 4
    #header = True
    #sep = '\t'
    #blast_identity = 0.7
    #cluster_identity = 0.9
    #dbscan_eps = 0.75
    
    #if not os.path.exists(out_dir):
    #    os.makedirs(out_dir)

    contigs = construct_sequences(fasta)
    read_xprs(contigs, xprs, n_col, x_col, header, '\t')
    
    print('    - start analyzing similar subsequences in contigs...')
    c_self_blastn = read_blastn(blastn, identity_t=blast_identity)
    c_ss_matches = find_match(c_self_blastn, contigs, copy.deepcopy(contigs))
    
    print('    - start building connected component for each sequence...')
    uf, clusters = connected_component(contigs, c_ss_matches, threshold=cluster_identity)
    
    print('    - start clustering xprs for each connected component...')
    quant_eval(clusters, dbscan_eps, out_file)
    
    print('finished, time elapsed: ', np.around(time.time() - start_time, 3), ' seconds', sep='')
    
if __name__ == '__main__':
    if len(sys.argv) != 12:
        print('usage: ', os.path.basename(__file__), 
              ' [blastn] [fasta] [xprs] [n_col] [x_col] [header] [sep] [blast identity] [cluster identity] [dbscan eps] [out file]', sep='')
        sys.exit(1)
    blastn = sys.argv[1]
    fasta = sys.argv[2]
    xprs = sys.argv[3]
    n_col = int(sys.argv[4])
    x_col = int(sys.argv[5])
    header = sys.argv[6]
    sep = sys.argv[7]
    blast_identity = float(sys.argv[8])
    cluster_identity = float(sys.argv[9])
    dbscan_eps = float(sys.argv[10])
    out_file = sys.argv[11]
    
    main(blastn, fasta, xprs, n_col, x_col, header, sep, blast_identity, cluster_identity, dbscan_eps, out_file)