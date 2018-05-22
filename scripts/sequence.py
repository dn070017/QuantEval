
from collections import defaultdict

class Sequence:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        
        self.label = dict()
        self.tr = dict()
        self.xprs = dict()
        self.contribute_xprs = defaultdict(dict)
        self.relative_xprs = defaultdict(dict)
            
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, target):
        return isinstance(target, Sequence) and self.name == target.name