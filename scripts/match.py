
import numpy as np

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
        
        identity = data['identity']
        
        alignment = np.zeros(q_seq.length).astype(float)
        alignment[q_start-1:q_end] = identity
        self.q_depth = alignment
        
        alignment = np.zeros(r_seq.length).astype(float)
        alignment[r_start-1:r_end] = identity
        self.r_depth = alignment
        
        self.q_global_identity = 0
        self.r_global_identity = 0
        
    def __hash__(self):
        return hash(self.m_name)
        
    def __eq__(self, target):
        return isinstance(target, Match) and self.m_name == target.m_name
    
    def extend(self, data, q_seq, r_seq):
        
        q_start = data['q_start']
        q_end = data['q_end']
        r_start = data['r_start']
        r_end = data['r_end']
        identity = data['identity']
        
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