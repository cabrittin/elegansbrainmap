"""
subcellular.py

Submodule for analyzing subcellular specificity

Author: Christopher Brittin
Created: 07 February 2018

"""


import networkx as nx
from random import shuffle

import aux
import db

class Subcell:
    def __init__(self,neuron,contin,endpts):
        self.neuron = neuron
        self.contin = contin
        self.endpts = endpts

        self.pre = {}
        self.post = {}
        self.gap = {}


    def generate_graph(self,cur,end=500):
        self.G = DB.mine.display2(cur,self.neuron)
        nodes = self.G.nodes()
        for n in nodes:
            if self.G.node[n]['contin'] != self.contin or self.G.node[n]['z'] > end:
                self.G.remove_node(n)

        self.length = nx.shortest_path_length(self.G,
                                              source=self.endpts[0],
                                              target=self.endpts[1])

        self.length = float(self.length)

    def generate_synapse_lists(self,cur):
        self.get_presynapses(cur)
        self.get_postsynapses(cur)
        self.get_gap_junctions(cur)
        
    def get_presynapses(self,cur):
        self.pre = self.get_synapses(cur,'post','preobj','chemical')
        self.pre_norm = self.compute_norm_dist(self.pre)

    def get_postsynapses(self,cur):
        post1 = self.get_synapses(cur,'pre','postobj1','chemical')
        post2 = self.get_synapses(cur,'pre','postobj2','chemical')
        post3 = self.get_synapses(cur,'pre','postobj3','chemical')
        post4 = self.get_synapses(cur,'pre','postobj4','chemical')

        self.post = self.add_dict(self.post,post1)
        self.post = self.add_dict(self.post,post2)
        self.post = self.add_dict(self.post,post3)
        self.post = self.add_dict(self.post,post4)
        self.post_norm = self.compute_norm_dist(self.post)

    def get_gap_junctions(self,cur):
        gap1 = self.get_synapses(cur,'post','preobj','electrical')
        gap2 = self.get_synapses(cur,'pre','postobj1','electrical')
        
        self.gap = self.add_dict(self.gap,gap1)
        self.gap = self.add_dict(self.gap,gap2)
        self.gap_norm = self.compute_norm_dist(self.gap)
        
        
    def get_synapses(self,cur,cellcol,objcol,stype):
        sql = ("select %s,%s,sections "
               "from synapsecombined "
               "join object "
               "on synapsecombined.%s = object.OBJ_Name "
               "where object.CON_Number = %d "
               "and synapsecombined.type = '%s'"
               %(objcol,cellcol,objcol,self.contin,stype))
        cur.execute(sql)

        tmp = {}
        for (obj,cell,weight) in cur.fetchall():
            #print(obj,cell,weight)
            if not self.G.has_node(int(obj)): continue
            cell = cell.split(',')
            for c in cell:
                if c not in tmp: tmp[c] = {}
                tmp[c][int(obj)] = int(weight)
        return tmp

    def random_synapses(self):
        shuffle(self.pre_obj)
        shuffle(self.post_obj)

        pre = self._random_synapses(self.pre,self.pre_obj)
        post = self._random_synapses(self.post,self.post_obj)

        pre = self.compute_norm_dist(pre)
        post = self.compute_norm_dist(post)

        return [pre,post]
        
    def _random_synapses(self,data,objs):
        idx = 0
        tmp = {}
        for n in data:
            tmp[n] = {}
            for o in data[n]:
                tmp[n][objs[idx]] = data[n][o]
                idx += 1
        return tmp
 
    def compute_norm_dist(self,data):
        dist = {}
        for n in data:
            num,den = 0.,0
            for o in data[n]:
                den += data[n][o]
                d = nx.shortest_path_length(self.G,source=self.endpts[0],target=o)
                num += data[n][o]*d
            if den > 0: dist[n] = num/den/self.length
        #return sorted(dist, key=lambda k: dist[k])
        return dist

    def _compute_dist(self,o):
        d = nx.shortest_path_length(self.G,source=self.endpts[0],target=o)
        return float(d)/self.length
    
    def add_dict(self,d1,d2):
        for n in d2:
            if n not in d1: d1[n] = {}
            for o in d2[n]:
                if o not in d1[n]: d1[n][o] = 0
                d1[n][o] += d2[n][o]
        return d1

    def initiate_randomize(self):
        self.pre_obj = self._initiate_randomize(self.pre)
        self.post_obj = self._initiate_randomize(self.post)
        
    def _initiate_randomize(self,data):
        objs = []
        for n in data:
            for o in data[n]:
                objs.append(o)
        return objs



