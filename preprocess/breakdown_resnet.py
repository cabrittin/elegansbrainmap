"""
@name: breakdown_resnet.py
@description:
    break down resnet synapses

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-04
"""
import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
import networkx as nx
from collections import defaultdict

from connectome.format_graphs import *
import ioaux

#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

def clean_graph(G,A):
    edges = [(a,b) for (a,b) in G.edges()]
    for (a,b) in edges:
        if not A.has_edge(a,b):
            G.remove_edge(a,b)

def get_ctype(l1,l2,b1,b2):
    ctype = -1
    if l2 > l1:
        ctype = ffwd_type(l1,l2,b1,b2)
    elif l2 < l1:
        ctype = fbk_type(l1,l2,b1,b2)
    else:
        ctype = rec_type(l1,l2,b1,b2)
    return ctype

def rec_type(l1,l2,b1,b2):
    ctype = -1
    if l1 == 1:
        ctype = 'sensory_inter_rec'
        if b1 == b2: ctype = 'sensory_rec'
    elif l1 == 2:
        ctype = 'l2_inter_rec'
        if b1 == b2: ctype = 'l2_rec'
    elif l1 == 3:
        ctype = 'l3_inter_rec'
        if b1 == b2: ctype = 'l3_rec'
    return ctype

def ffwd_type(l1,l2,b1,b2):
    ctype = -1
    if b1 == b2:
        ctype = 'ffwd_norm'
        if l2 - l1 > 1: ctype = 'ffwd_skip'
    else:
        ctype = 'ffwd_violation'
        if l1 == 1 and l2 == 3: ctype='cross_sensory_skip'
    return ctype

def fbk_type(l1,l2,b1,b2):
    ctype = -1
    if b1 == b2:
        ctype = 'fbk'
        if l1 - l2 > 1: ctype = 'fbk_skip'
    else:
        ctype = 'fbk_violation'
    return ctype

rm = ['PVR','HSNL','HSNR','PLNL','PLNR','PVNL','PVNR','PVR.']
RESNET = ['ffwd_norm','ffwd_skip','l2_rec','l3_inter_rec','l3_rec','sensory_rec']
FBK = ['fbk','fbk_skip']


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('-g','--group',
                        action = 'store',
                        dest = 'group',
                        required = False,
                        default = None,
                        help = 'Axon group'
                        )
   
    params = parser.parse_args()

    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    
    bundles = ioaux.read.into_dict(cfg['clusters']['brainmap_lr'])
    resnet = ioaux.read.into_dict(cfg['mat']['resnet']) 
    nbl = dict([(k,(v,int(resnet[k]))) for (k,v) in bundles.items() if int(resnet[k]) > -1])

    A = nx.read_graphml(cfg['refgraphs']['adj']%4)
    C = nx.read_graphml(cfg['refgraphs']['chem']%4)
    A.remove_nodes_from(rm)
    C.remove_nodes_from(rm)
    
    clean_graph(C,A)
    
    SD = [['pre','post','in_resnet','in_fbk']]

    ctype = defaultdict(int)
    cweight = defaultdict(int)
    resedges = []
    for (u,v) in C.edges():
        (b1,l1) = nbl[u]
        (b2,l2) = nbl[v]
        w = C[u][v]['weight']
        _ctype = get_ctype(l1,l2,b1,b2)
        resedges.append([u,v,_ctype])
        in_resnet = int(_ctype in RESNET)
        in_fbk = int(_ctype in FBK)
        SD.append([u,v,in_resnet,in_fbk])
        if _ctype == params.group: print(u,v,w)
        ctype[_ctype] += 1
        cweight[_ctype] += w

    print('Aggregate of synaptic contacts')
    total = 0
    for (k,v) in sorted(ctype.items()):
        print('%s,%d'%(k,v))
        total += v
    
    print('\n\n')
    print('Aggregate of synaptic EM sections')
    for (k,v) in sorted(cweight.items()):
        print('%s,%d'%(k,v))
        total += v


    print(total)
    ioaux.write.from_list(cfg['resnet']['resnet_edges'],resedges)
    ioaux.write.from_list('source_data/source_data_figure_4a.csv',SD)


    
