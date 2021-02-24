"""
@name: prune_polyads.py
@description:
Make prune C1 polyads from reference graphs.
For a given polyadic synapse the C1 postsynaptic partners if there is at least one C2 or greater contact.

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-05
"""

import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import networkx as nx

from connectome.format_graphs import make_reference_graphs,clean_graph
from connectome.load import from_db,filter_synapses
from connectome.format_graphs import consensus_graph,filter_graph_edge,normalize_edge_weight
from connectome.load import reference_graphs
import aux


#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

def sfilter(synapses,args):
    Ref,lrmap,left,right = args
    pruned = []
    count = 0
    for s in synapses:
        reflected = False
        pre = s[0]
        post = s[1].split(',')
        if pre in right:
            reflected = True            
            pre = lrmap[pre]
            _post = []
            for p in post:
                if p in lrmap: _post.append(lrmap[p])
            post = _post
        c1 = []
        c2 = []
        for p in post:
            if Ref.has_edge(pre,p):
                if Ref[pre][p]['id'] == 1:
                    c1.append(p)
                else:
                    c2.append(p)
        if c1 and c2:
            post = c2
            if reflected: post = [lrmap[p] for p in post]
            s = (s[0],','.join(post),s[2],s[3],s[4])
            count += len(c1)
        pruned.append(s)
    #print(count)
    return pruned

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    

    params = parser.parse_args()

    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)

    left = aux.read.into_list(cfg['mat']['left_nodes'])
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = aux.read.into_lr_dict(cfg['mat']['lrmap'])
    nodes = aux.read.into_list(cfg['mat']['nodes']) 
    remove = aux.read.into_list(cfg['mat']['remove'])
    edge_thresh = cfg.getint('params','lower_weight_threshold')
    dbs = cfg['input']['databases'].split(',')

    A4 = nx.read_graphml(cfg['refgraphs']['adj']%4)
    C = dict([(i,nx.read_graphml(cfg['refgraphs']['chem']%i)) for i in range(1,5)])
    Ref = make_reference_graphs(C)
    Ref = clean_graph(Ref,A4)

    G = []
    for d in ['JSH','N2U']:#dbs:
        D = from_db(d,adjacency=True,chemical=True,electrical=True,remove=remove,dataType='networkx')
        P = filter_synapses(d,sfilter,args=[Ref,lrmap,left,right],remove=remove)
        P.C = P.C.subgraph(D.C.nodes())
        print(D.C.number_of_edges(),P.C.number_of_edges())
        D.C = P.C
        D.A = filter_graph_edge(D.A,pct=edge_thresh)
        D.split_left_right(left,right)  
        D.map_right_graphs(lrmap)
        G.append(D)
    
    total,_total = 0,0
    for i in range(1,5):
        pC = nx.DiGraph()
        #consensus_graph(pC,[G[0].Cl,G[0].Cr,G[1].Cl,G[1].Cr],i,left)
        consensus_graph(pC,[G[0].Cl,G[0].Cr,G[1].Cl,G[1].Cr],i,nodes)
        print(pC.number_of_edges(),C[i].number_of_edges())
        pC = clean_graph(pC,A4)
        cC = clean_graph(C[i],A4)
        print(i,pC.number_of_edges(),cC.number_of_edges())
        total += pC.number_of_edges()
        _total += cC.number_of_edges()
        nx.write_graphml(pC,cfg['refgraphs']['pruned_chem']%i)

    print(total,_total)
