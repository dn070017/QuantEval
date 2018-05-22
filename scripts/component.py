
import numpy as np

from collections import defaultdict

class Component:
    def __init__(self, name):
        self.name = name
        self.member = list()
        self.member_xprs = defaultdict(list)

    def add_member(self, seq):
        self.member.append(seq)
        
        for xprs_label in seq.xprs:
            self.member_xprs[xprs_label].append(seq.xprs[xprs_label])
        
    def get_maximum_xprs(self, label):
        return round(np.amax(self.member_xprs[label]), 3)
    
    def get_total_xprs(self, label):
        return round(np.sum(self.member_xprs[label]), 3)
    
    def get_average_xprs(self, label):
        return round(np.average(self.member_xprs[label]), 3)
    
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, target):
        return isinstance(target, Component) and self.name == target.name
    
class Cluster(Component):
    def __init__(self, name, xprs_label):
        self.name = name
        self.xprs_label = xprs_label
        
        self.member = list()
        self.member_xprs = list()
        
        self.rank = np.nan
        
    def add_member(self, seq):
        self.member.append(seq)
        self.member_xprs.append(seq.xprs[self.xprs_label])
        
    def get_maximum_xprs(self):
        return round(np.amax(self.member_xprs), 3)
    
    def get_total_xprs(self):
        return round(np.sum(self.member_xprs), 3)
    
    def get_average_xprs(self):
        return round(np.average(self.member_xprs), 3)
    
    def __gt__(self, target):
        return isinstance(target, Cluster) and self.get_average_xprs() <= target.get_average_xprs()