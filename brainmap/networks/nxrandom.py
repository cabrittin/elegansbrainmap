"""
@name:nxrandom.py

Class for generated randomized networks
Assumes networkx input

@description

@author: Christopher Brittin
@email: 'cabrittin' + <at> + 'gmail' + '.' + 'com'
"""

import argparse
import networkx as nx
import random
import numpy as np

class Randomize(nx.Graph):
    def __init__(self,G):
        if G.is_directed():
            self.R = nx.DiGraph()
        else:
            self.R = nx.Graph()
        for (u,v,w) in G.edges.data(data='weight'): self.R.add_edge(u,v,weight=w)
        self.num_switches = 0

    def edge_switch(self,iters=1000):
        self.num_switches = 0
        while self.num_switches < iters:
            [(u1,v1),(u2,v2)] = random.sample(self.R.edges,2)
            cond1 = len(set([u1,v1,u2,v2])) == 4
            cond2 = not self.R.has_edge(u1,v2)
            cond3 = not self.R.has_edge(u2,v1)
            cond4 = self.check_cond4((u1,v1),(u2,v2))
            if cond1 and cond2 and cond3 and cond4:
                w1 = self.R[u1][v1]['weight']
                w2 = self.R[u2][v2]['weight']
                self.R.remove_edges_from([(u1,v1),(u2,v2)])
                self.R.add_edge(u1,v2,weight=w1)
                self.R.add_edge(u2,v1,weight=w2)
                self.num_switches += 1
    
    def check_degree_dist(self,A):
        disc = 0
        disc_nodes = []
        for n in A.nodes():
            if A.degree(n) != self.R.degree(n):
                #print(n,A.degree(n),self.degree(n))
                disc += 1
                disc_nodes.append(n)
        return disc,disc_nodes 

    def check_similarity(self,A):
        n = A.number_of_nodes()
        sim = np.zeros(n)
        for (i,m) in enumerate(A.nodes()):
            neigh1 = self.R.neighbors(m)
            neigh2 = A.neighbors(m)
            sim[i] = self.similarity_measure(neigh1,neigh2)
        return np.mean(sim)

    def similarity_measure(self,neigh1,neigh2):
        neigh1 = set(neigh1)
        neigh2 = set(neigh2)
        num = float(len(neigh1 & neigh2))
        den = len(neigh1 | neigh2)
        if den: 
            return num / den
        else:
            return 0
   
    def load_condition_params(self):
        pass

    def check_cond4(self,e1,e2):
        return True
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fin',
                        action = 'store',
                        help = 'Input graphml file')
    params = parser.parse_args()

    A = nx.read_graphml(params.fin)

    print(A.number_of_edges())

    R = Randomize(A)
    print(A.number_of_nodes(),R.R.number_of_nodes())
    print(A.number_of_edges(),R.R.number_of_edges())
    #R.edge_switch(iters=1000)
    print('Number of switches: %d' %R.num_switches)
    disc,disc_nodes = R.check_degree_dist(A)
    print('Total number of discrepant degrees: %d' %disc)
    sim =  R.check_similarity(A)
    print('Mean similarity across nodes: %1.3f' %sim)
