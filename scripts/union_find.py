
class UnionFind:
   
    def __init__(self, seq_dict):
        self.parent = dict()
        self.size = dict()
        self.component_label = dict()
        self.component_size = dict()
        
        for name in seq_dict.keys():
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
    
    def get_component_size(self, c):
        return self.component_size[c]
    
    def rename_component(self):
        index = 1
        for name, parent in self.parent.items():
            if name == parent:
                label = 'component_' + "{0:010d}".format(index)
                self.component_label[name] = label
                self.component_size[label] = self.size[name]
                index += 1
        
        for name in self.parent.keys():
            root = self.find(name)
            self.component_label[name] = self.component_label[root]
    
        return