"""
@name: model_init.py

@description
Format data for  3-param model.


@author: Christopher Brittin
@email: 'cabrittin' + <at> + 'gmail' + '.' + 'com'
"""

import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation   
import networkx as nx
import numpy as np

from connectome.format_graphs import *
from connectome.load import reference_graphs
import aux

def add_intra_inter_bundle(Ref,bundles):
    for (a,b) in Ref.edges():
        Ref[a][b]['intra'] = 1
        try:
            if bundles[a] != bundles[b]: Ref[a][b]['intra'] = 0
        except:
            pass

def breakdown_count(Ref):
    count = np.zeros((3,4))
    for (a,b,v) in Ref.edges.data(data=True):
        idx = v['id'] - 1
        count[0,idx] += 1
        if v['intra']:
            count[1,idx] += 1
        else:
            count[2,idx] += 1 
    return count

def high_pass_standardized_filter(G,thresh):
    edges = [(a,b) for (a,b,w) in G.edges.data('weight') if w < thresh]
    G.remove_edges_from(edges)

def low_pass_standardized_filter(G,thresh):
    edges = [(a,b) for (a,b,w) in G.edges.data('weight') if w >= thresh]
    G.remove_edges_from(edges)

def head_edges(Ref,limit=10):
        counter = 0
        for e in Ref.edges.data(): 
            if counter == limit: break
            print(i,e)
            counter += 1
 

CONFIG = os.environ['CONFIG']


def run(params):
    BAD_NODES = ['PVR','PVNL','HSNL','PVR.','HSNR','SABD','VB02','PLNL','PLNR','PVNR']
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    
    #bundles = aux.read.into_dict(cfg['clusters']['m4'])
    bundles = aux.read.into_dict(cfg['clusters']['final'])
    #intra_bundle = cfg.getboolean('model_params','intra_bundle')
    #weight_screen = cfg.getboolean('model_params','weight_screen')
    restricted=cfg.getboolean('model_params','restricted') 
    clean_idx = params.delta
    
    if not restricted: BAD_NODES = []

    A,C,E = reference_graphs(cfg)
    
    print('#M4',A[4].number_of_nodes())
    Ref = make_reference_graphs(A)
    lRef = Ref.copy()
    hRef = Ref.copy()
    standardize_edge_weigth(lRef)
    standardize_edge_weigth(hRef)
    low_pass_standardized_filter(lRef,0)
    high_pass_standardized_filter(hRef,0)

    for i in range(1,5): print(A[i].number_of_edges())
    data = [A,C,E]
    all_count = np.zeros((3,4))
    low_count = np.zeros((3,4))
    high_count = np.zeros((3,4))
    intra_count = np.zeros((3,4))
    inter_count = np.zeros((3,4))

    for (i,G) in enumerate(data):  
        Ref = make_reference_graphs(G)
        #head_edges(Ref) 
        for n in BAD_NODES:
            if Ref.has_node(n): Ref.remove_node(n)
        if i > 0: Ref = clean_graph(Ref,A[clean_idx]) #only clean Chem and Gap graphs
        ###If removing 1 EM sections -- othwerwise comment out###
        #if i == 1:
        #    rm_edge = [(a,b) for (a,b) in Ref.edges() if Ref[a][b]['weight'] < 2]
        #    Ref.remove_edges_from(rm_edge)
        #head_edges(Ref) 
        add_intra_inter_bundle(Ref,bundles)
        count = breakdown_count(Ref) 
        all_count[i,:] = count[0,:]
        intra_count[i,:] = count[1,:]
        inter_count[i,:] = count[2,:]
        
        #Weight screen
        _lRef = Ref.copy()
        _hRef = Ref.copy()
        #head_edges(_lRef)
        print('1',_lRef.number_of_edges())
        _lRef.remove_edges_from([e for e in _lRef.edges() if not lRef.has_edge(e[0],e[1])])
        _hRef.remove_edges_from([e for e in _hRef.edges() if not hRef.has_edge(e[0],e[1])])
        print('2',_lRef.number_of_edges())
        #standardize_edge_weigth(lRef)
        #standardize_edge_weigth(hRef)
        #low_pass_standardized_filter(lRef,0)
        #high_pass_standardized_filter(hRef,0)
        #add_intra_inter_bundle(lRef,bundles)
        #add_intra_inter_bundle(hRef,bundles)
        lcount = breakdown_count(_lRef)
        hcount = breakdown_count(_hRef)
        low_count[i,:] = lcount[0,:]
        high_count[i,:] = hcount[0,:]

    print(all_count)
    print(low_count)
    print(high_count)
    print(A[4].number_of_nodes())
    print(all_count / all_count.sum(axis=1,keepdims=True))
       
    np.save(cfg['model_data']['all'],all_count)
    #np.save(cfg['model_data']['all'].replace('.npy','_%d.npy'%params.delta),all_count)
    np.save(cfg['model_data']['inter_bundle'],inter_count)
    np.save(cfg['model_data']['intra_bundle'],intra_count)
    np.save(cfg['model_data']['low_weight'],low_count)
    np.save(cfg['model_data']['high_weight'],high_count)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('-d','--delta',
                        dest = 'delta',
                        action = 'store',
                        default = 4,
                        type = int,
                        required = False,
                        help = 'Delta value')
 
    params = parser.parse_args()

    run(params)
