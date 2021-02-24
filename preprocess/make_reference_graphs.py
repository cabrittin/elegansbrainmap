"""
@name: make_reference_graphs.py
@description:
Make reference graphs
Modified for rosettes


@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-05
"""

import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import networkx as nx
import time

from connectome.load import from_db
from connectome.format_graphs import consensus_graph,filter_graph_edge,normalize_edge_weight
import aux


#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini' 

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('deg',
                        action = 'store',
                        type = int,
                        help = 'Degree of reproducibility')
 
    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('--adj',
                        dest = 'adj',
                        action = 'store_true',
                        default = False,
                        required = False,
                        help = 'Display axon contacts')
     
    parser.add_argument('--chem',
                        dest = 'chem',
                        action = 'store_true',
                        default = False,
                        required = False,
                        help = 'Display chemical contacts')
    
    parser.add_argument('--gap',
                        dest = 'gap',
                        action = 'store_true',
                        default = False,
                        required = False,
                        help = 'Display gap junction contacts')


    params = parser.parse_args()

    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    deg = params.deg

    left = aux.read.into_list(cfg['mat']['left_nodes'])
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = aux.read.into_lr_dict(cfg['mat']['lrmap'])
    nodes = aux.read.into_list(cfg['mat']['nodes'])
    remove = aux.read.into_list(cfg['mat']['remove'])
    edge_thresh = cfg.getint('params','lower_weight_threshold')
    dbs = cfg['input']['databases'].split(',')
    
    G = []
    for d in dbs:
        D = from_db(d,adjacency=True,chemical=True,electrical=True,remove=remove,dataType='networkx')
        D.A = filter_graph_edge(D.A,pct=edge_thresh)
        D.split_left_right(left,right)  
        D.map_right_graphs(lrmap)
        G.append(D)
   

    if params.adj:
        A = nx.Graph()
        for g in G: 
            normalize_edge_weight(g.Al)
            normalize_edge_weight(g.Ar)
        #consensus_graph_deprecated(A,[G[0].Al,G[0].Ar,G[1].Al,G[1].Ar],deg,left)
        consensus_graph(A,[G[0].Al,G[0].Ar,G[1].Al,G[1].Ar],deg,nodes,weight=['weight','wnorm'])
        #for e in A.edges.data(): print(e)
        nx.write_graphml(A,cfg['refgraphs']['adj']%deg)
        print(f"Wrote to : {cfg['refgraphs']['adj']%deg}")
    
    if params.chem: 
        C = nx.DiGraph()
        #consensus_graph(C,[G[0].Cl,G[0].Cr,G[1].Cl,G[1].Cr],deg,left)
        consensus_graph(C,[G[0].Cl,G[0].Cr,G[1].Cl,G[1].Cr],deg,nodes)
        #for e in C.edges.data(): print(e)
        nx.write_graphml(C,cfg['refgraphs']['chem']%deg)
        print(f"Wrote to : {cfg['refgraphs']['chem']%deg}")
    
    if params.gap: 
        E = nx.Graph()
        #consensus_graph(E,[G[0].El,G[0].Er,G[1].El,G[1].Er],deg,left)
        consensus_graph(E,[G[0].El,G[0].Er,G[1].El,G[1].Er],deg,nodes)
        nx.write_graphml(E,cfg['refgraphs']['gap']%deg)
        print(f"Wrote to : {cfg['refgraphs']['chem']%deg}")
