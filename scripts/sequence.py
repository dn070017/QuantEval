
import numpy as np

class Sequence:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        
        self.label = dict()
        self.size = dict()
        self.tr = dict()
        self.xprs = dict()
        
        #self.tr_good = np.nan
        #self.tr_bases_covered = np.nan 
        #self.tr_seq_true = np.nan 
        #self.tr_score =np.nan 
        #self.tr_not_segmented = np.nan        
        
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, target):
        return isinstance(target, Sequence) and self.name == target.name