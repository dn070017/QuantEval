#!/usr/bin/env python3

import numpy as np
import os
import pickle
import pandas as pd
import sys

def usage():
    terminate = False
    
    if len(sys.argv) != 6:
        terminate = True
    else:
        if sys.argv[1] not in ['simulation', 'real']:
            terminate = True
    if terminate:
        print('usage:', os.path.basename(__file__), '[simulateion/real] [flux simulator pro/kallisto]' +
              ' [flux simulator lib/rsem] [unpaired pickle/salmon] [output tsv]')
        sys.exit(1)
    
    return

def eff_len(L, pmf):
    min_val = 0
    max_val = len(pmf)-1
    c_len = min_val
    max_length = min(L, max_val)
    effective_length = 0
    while (c_len <= max_length):
        i = c_len - min_val
        effective_length += pmf[i] * (L - c_len + 1)
        c_len += 1
    return effective_length

def read_lib(lib):
    count_dict = dict()
    with open(lib) as lib_in:
        for l in lib_in:
            data = l.rstrip().split()
            start, end = int(data[0]), int(data[1])
            length = abs(start - end)
            count = int(data[-1])
            if length in count_dict:
                count_dict[length] += count
            else:
                count_dict[length] = count
    return count_dict

def estimate_eff_length(pro, lib):
    count_dict = read_lib(lib)
    length_list = list()
    count_list = list()
    
    for length, count in count_dict.items():
        length_list.append(length)
        count_list.append(count)
    total_count = float(sum(count_list))
    max_frag = max(length_list)
    pdf = np.zeros(max_frag + 1)
    cdf = np.zeros(max_frag + 1)
    for length, count in zip(length_list, count_list):
        pdf[length] += count / total_count
    cdf = np.cumsum(pdf)
    
    vfunc = np.vectorize(eff_len, excluded=['pmf'])
    mean_length = pdf.mean()
    # Assuming the current lengths are in d["Length_true"]
    original_length = pro["length"].values
    eff_length = vfunc(L=original_length, pmf=pdf)
    pro["eff_length"] = eff_length
    
def calculate_answer_tpm(flux_simulator_pro, flux_simulator_lib, unpaired_pickle, output):
    flux_columns = ['locus', 'name', 'conding', 'length', 'molecular_fraction', 'molecular_count', 
                    'fragment_fraction', 'fragment_count', 'read_fraction', 'read_count', 'covered', 
                    'chi_square', 'variation']
    flux_simulator = pd.read_table(flux_simulator_pro, sep='\t', header=None)
    flux_simulator.columns = flux_columns
    
    estimate_eff_length(flux_simulator, flux_simulator_lib)

    with open(unpaired_pickle, 'rb') as count_pickle:
        unpaired_read_count = pickle.load(count_pickle)
        
    unpaired_read_count = pd.DataFrame.from_dict(unpaired_read_count, orient='index')
    unpaired_read_count.columns = ['count']
    unpaired_read_count['name'] = unpaired_read_count.index.tolist()

    flux_simulator = pd.merge(unpaired_read_count, flux_simulator, on='name', how='outer') 
    flux_simulator = flux_simulator.fillna(value=0)
    flux_simulator['read_per_nucleotide'] = (flux_simulator['read_count'] - flux_simulator['count']) / flux_simulator['eff_length']
    flux_simulator['answer_tpm'] = 10 ** 6 * flux_simulator['read_per_nucleotide'] / np.sum(flux_simulator['read_per_nucleotide'])
    flux_simulator['answer_count'] = flux_simulator['read_count']
    flux_simulator = flux_simulator.loc[:, ('name', 'answer_tpm', 'answer_count')]
    flux_simulator = flux_simulator.round({'answer_tpm': 3})

    flux_simulator.to_csv(output, sep='\t', index=False)
    
def average_tpm(kallisto, rsem, salmon, output):
    kallisto = pd.read_table(kallisto, sep='\t')
    kallisto = kallisto.rename(columns={'target_id': 'name', 'tpm': 'kallisto_tpm', 'est_counts': 'kallisto_count'})
    kallisto = kallisto.loc[:, ('name', 'kallisto_tpm', 'kallisto_count')]

    rsem = pd.read_table(rsem, sep='\t')
    rsem = rsem.rename(columns={'transcript_id': 'name', 'TPM': 'rsem_tpm', 'expected_count': 'rsem_count'})
    rsem = rsem.loc[:, ('name', 'rsem_tpm', 'rsem_count')]

    salmon = pd.read_table(salmon, sep='\t')
    salmon = salmon.rename(columns={'Name': 'name', 'TPM': 'salmon_tpm', 'NumReads': 'salmon_count'})
    salmon = salmon.loc[:, ('name', 'salmon_tpm', 'salmon_count')]

    tmp = pd.merge(salmon, rsem, on='name', how='inner')
    tmp = pd.merge(kallisto, tmp, on='name', how='inner')
    tmp['answer_tpm'] = (tmp['kallisto_tpm'] + tmp['rsem_tpm'] + tmp['salmon_tpm']) / 3
    tmp['answer_count'] = (tmp['kallisto_count'] + tmp['rsem_count'] + tmp['salmon_count']) / 3
    expression_table = tmp.loc[:, ('name', 'answer_tpm', 'answer_count')]
    expression_table = expression_table.round({'answer_tpm': 3})

    expression_table.to_csv(output, sep='\t', index=False)
    
if __name__ == '__main__':
    usage()
    if sys.argv[1] == 'simulation':
        calculate_answer_tpm(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        average_tpm(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])