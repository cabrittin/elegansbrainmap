import numpy as np
from scipy.special import digamma

def bin_for_knn(w):
    if w > 30:
        return 10
    elif w > 20:
        return 9
    elif w > 15:
        return 8
    elif w > 10:
        return 7
    elif w > 7:
        return 6
    elif w > 5:
        return 5
    else:
        return int(w-1)
    
def get_groups(data):
    groups = {}
    for i in set(data[:,1]):
        groups[i] = Group(i)
        jdx = np.where(data[:,1]==i)[0]
        for j in jdx:
            groups[i].add_vertex(Vertex(j))
        if groups[i].Nx > 1: groups[i].vertices[0].add_right(groups[i].vertices[1])
        if groups[i].Nx > 1: groups[i].vertices[-1].add_left(groups[i].vertices[-2])
        for j in range(1,groups[i].Nx - 1):
            groups[i].vertices[j].add_left(groups[i].vertices[j-1])
            groups[i].vertices[j].add_right(groups[i].vertices[j+1])
            
    return groups

def compute_mutual_info(groups,kradius):
    psiNx = []
    psim = []
    N = 0
    for i in groups:
        N += groups[i].Nx
        for v in groups[i].vertices:
            v.knn_search(kradius)
            psiNx.append(digamma(groups[i].Nx))
            psim.append(digamma(v.m))
            #print(i,v.id,v.bounds,v.m)

    psiNx = np.mean(psiNx)
    psim = np.mean(psim)
    psik = digamma(kradius)
    psiN = digamma(N)
    I = psiN - psiNx + psik - psim    

    return I
    
class Group:
    def __init__(self,group):
        self.group = group
        self.Nx = 0
        self.vertices = []

    def add_vertex(self,id):
        self.Nx += 1
        self.vertices.append(id)

class Vertex:
    def __init__(self,id):
        self.id = id
        self.left_v = None
        self.left_d = 0
        self.right_v = None
        self.right_d = 0
        self.bounds = [None,None]

    def add_left(self,vertex):
        self.left_v = vertex
        self.left_d = abs(vertex.id - self.id)

    def add_right(self,vertex):
        self.right_v = vertex
        self.right_d = abs(vertex.id - self.id)

    def get_left(self):
        return self.left_v,self.left_d

    def get_right(self):
        return self.right_v,self.right_d

    def knn_search(self,kradius):
        if kradius == 0: return []
        dist = []
        _d = 0
        _kradius = kradius
        v,d = self.get_left()
        while v and _kradius:
            _d += d
            dist.append((v.id,_d))
            v,d = v.get_left()
            _kradius -= 1

        _d = 0
        _kradius = kradius
        v,d = self.get_right()
        while v and _kradius:
            _d += d
            dist.append((v.id,_d))
            v,d = v.get_right()
            _kradius -= 1

        dist = sorted(dist, key=lambda x: x[1])[:kradius]
        _bounds = sorted([self.id] + [a[0] for a in dist])
        self.bounds = [_bounds[0],_bounds[-1]]
        self.m = self.bounds[1] - self.bounds[0] + 1
