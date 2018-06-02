
from collections import defaultdict

class Sequence:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        
        self.label = dict()
        self.tr = dict()
        self.xprs_tpm = dict()
        self.relative_xprs_tpm = defaultdict(dict)
        self.contribute_xprs_tpm = defaultdict(dict)
        
        self.xprs_count = dict()
        self.relative_xprs_count = defaultdict(dict)
        self.contribute_xprs_count = defaultdict(dict)
            
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, target):
        return isinstance(target, Sequence) and self.name == target.name