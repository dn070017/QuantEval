
import numpy as np

from collections import defaultdict

class Component:
    def __init__(self, name):
        self.name = name
        self.size = 0
        self.member = list()
        self.xprs_tpm = defaultdict(list)
        self.xprs_count = defaultdict(list)

    def add_member(self, seq):
        self.member.append(seq)
        self.size += 1
        
        for xprs_label in seq.xprs_tpm:
            self.xprs_tpm[xprs_label].append(seq.xprs_tpm[xprs_label])
            
        for xprs_label in seq.xprs_count:
            self.xprs_count[xprs_label].append(seq.xprs_count[xprs_label])
        
    def get_maximum_xprs(self, metric, label):
        if metric == 'tpm':
            return round(np.amax(self.xprs_tpm[label]), 3)
        elif metric == 'count':
            return round(np.amax(self.xprs_count[label]), 3)
    
    def get_total_xprs(self, metric, label):
        if metric == 'tpm':
            return round(np.sum(self.xprs_tpm[label]), 3)
        elif metric == 'count':
            return round(np.sum(self.xprs_count[label]), 3)
    
    def get_average_xprs(self, metric, label):
        if metric == 'tpm':
            return round(np.average(self.xprs_tpm[label]), 3)
        elif metric == 'count':
            return round(np.average(self.xprs_count[label]), 3)
    
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, target):
        return isinstance(target, Component) and self.name == target.name