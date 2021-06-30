"""
@name: collapse_graph_nodes.py
@description:
Combines left/right nodes. By default, right nodes are 'absorbed' by the left nodes

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-05
"""
import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
import networkx as nx
import numpy as np

import ioaux

#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('fin',
                        action='store',
                        help = 'Input graphml file')

    parser.add_argument('--keep_self_loops',
                        action='store_true',
                        dest='self_loop',
                        required=False,
                        default=False,
                        help="Keep self loops, default is False")
    
    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
 
    params = parser.parse_args()
   
    fout = params.fin.replace('.graphml','_cl.graphml')
    
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    
    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    rlmap = dict(zip(right,left))

    G = nx.read_graphml(params.fin)
    if G.is_directed():
        H = nx.DiGraph()
    else:
        H = nx.Graph()

    nodes = set(G.nodes())
    lr = set(left + right)
    dif = nodes - lr
    
    for d in dif: rlmap[d] = d
    for l in left: rlmap[l] = l

    for (u,v) in G.edges:
        for (k,w) in G[u][v].items():
            um = rlmap[u]
            vm = rlmap[v]
            if not H.has_edge(um,vm): H.add_edge(um,vm)
            try:
                H[um][vm][k] += w
            except:
                H[um][vm][k] = w


    if not params.self_loop:
        H.remove_edges_from(nx.selfloop_edges(H,data=True))

    nx.write_graphml(H,fout)
    print(f"Collapse nodes from {G.number_of_nodes()} to {H.number_of_nodes()}")
    print(f"Writing to: {fout}")
