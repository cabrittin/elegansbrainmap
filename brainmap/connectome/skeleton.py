import numpy as np
import networkx as nx
from itertools import combinations



class Skeleton(nx.Graph):
    def __init__(self,contin):
        self.contin = contin
        nx.Graph.__init__(self)
        self.scale = np.array([5,5,18])

    def set_endpoints(self):
        epts = [v for v in self.nodes if self.degree(v) == 1]
        comb = combinations(epts,2)
        maxdist = 0
        maxepts = [0,0]
        for (a,b) in comb:
            try:
                length = nx.shortest_path_length(self,a,b,weight='weight')
                if length > maxdist:
                    maxdist = length
                    maxepts = [a,b]
            except:
                continue

        self.endpts = maxepts
        self.length = maxdist

    def set_distances(self):
        tmp = {}
        for (a,b) in self.edges():
            dx = np.multiply((self.node[a]['loc'] - self.node[b]['loc']),self.scale)
            dist = np.sqrt(np.sum(dx**2))
            self._adj[a][b]['weight'] = dist    

    def get_node_distance(self,target):
        dist = nx.shortest_path_length(self,self.endpts[0],target,weight='weight')
        return dist / self.length

    def plot_skeleton(self,ax):
        for (i,j) in self.edges():
            x = [self.node[i]['loc'][0],self.node[j]['loc'][0]]
            y = [self.node[i]['loc'][1],self.node[j]['loc'][1]]
            z = [self.node[i]['loc'][2]*self.scale[2],self.node[j]['loc'][2]*self.scale[2]]
            ax.plot(z,x,y,'b-')

    def plot_max_endpoints(self,ax):
        for i in self.endpts:
            x = self.node[i]['loc'][0]
            y = self.node[i]['loc'][1]
            z = self.node[i]['loc'][2]*self.scale[2]
            ax.plot([z],[x],[y],'ro')
            ax.text(z,x,y,str(i),(1,1,0))
        
    def plot_endpoints(self,ax):
        epts = [v for v in self.nodes if self.degree(v) == 1 and v not in self.endpts]
        for i in epts:
            x = self.node[i]['loc'][0]
            y = self.node[i]['loc'][1]
            z = self.node[i]['loc'][2]*self.scale[2]
            ax.plot([z],[x],[y],'go')
            ax.text(z,x,y,str(i),(1,1,0))
            


