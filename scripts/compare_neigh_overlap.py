"""
compare_neigh_overlap.py

Plots distributions of Jaccard distances for overlapping ipsilateral 
neighborhoods (blue) and homologous contralateral neighborhoods (red) 
in the adult and L4.

crated: Christopher Brittin
data: 01 November 2018

"""
import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Brittin modules
from connectome.load import from_db
from connectome.format_graphs import low_pass_edge_filter,mid_pass_edge_filter,high_pass_edge_filter
from networks.stats import get_neighborhood_similarity,get_neighborhood_overlap_similarity
import aux

CONFIG = os.environ['CONFIG']
SOURCE = "data/neighborhood_similarity.csv"

def add_neighborhood_similarity(data,label,A,reflected,left,right):
    for (a,b) in get_neighborhood_similarity(A,reflected,left):
        if b == -1 or b > 1: continue
        data.append(label + ['homologous',a,b])

def add_overlap_similarity(data,label,A,vertices):
    for (a,b) in get_neighborhood_overlap_similarity(A,vertices):
        if b == -1 or b > 1: continue
        data.append(label + ['proximal',a,b])

def ipsilateral_pass_filter(G,args):
    [left,right] = args
    edges = []
    for e in G.es:
        u = G.vs[e.source]['name']
        v = G.vs[e.target]['name']
        c1 = (u in left) and (v in right)
        c2 = (u in right) and (v in left)
        if c1 or c2: edges.append(e)
    G.delete_edges(edges)

def contralateral_pass_filter(G,args):
    [left,right] = args
    edges = []
    for e in G.es:
        u = G.vs[e.source]['name']
        v = G.vs[e.target]['name']
        c1 = (u in left) and (v in right)
        c2 = (u in right) and (v in left)
        if not (c1 or c2): edges.append(e)
    G.delete_edges(edges)

def add_data(cfg,data,_label,edge_filter=None,args=None):
    N2U = 'N2U'
    JSH = 'JSH'
    
    left = aux.read.into_list(cfg['mat']['left_nodes'])
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = aux.read.into_lr_dict(cfg['mat']['lrmap']) 
    #_remove = ['VC01','VD01','VB01','VB02','HSNL','HSNR','PVNL','PVNR']
    _remove = ['VC01','VD01','VB01','VB02','HSNL','HSNR','PVNL','PVNR','PLNL','PLNR','PVR','PVR.']

    label = ['Adult L/R'] + _label   
    n2u = from_db(N2U,adjacency=True,remove=_remove)
    if edge_filter: edge_filter(n2u.A,args)
    reflected = n2u.A.map_vertex_names(lrmap)
    add_neighborhood_similarity(data,label,n2u.A,reflected,left,right)
    add_overlap_similarity(data,label,n2u.A,left + right)
    
    label = ['L4 L/R'] + _label   
    jsh = from_db(JSH,adjacency=True,remove=_remove)
    if edge_filter: edge_filter(jsh.A,args)
    reflected = jsh.A.map_vertex_names(lrmap)
    add_neighborhood_similarity(data,label,jsh.A,reflected,left,right)
    add_overlap_similarity(data,label,jsh.A,left + right)
    
    label = ['Adult/L4'] + _label
    vertices = sorted((set(n2u.neurons)&set(jsh.neurons))-set(_remove))
    add_neighborhood_similarity(data,label,n2u.A,jsh.A,left,right)
    add_overlap_similarity(data,label,n2u.A,vertices)
    add_overlap_similarity(data,label,jsh.A,vertices)
  

def run(_cfg,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)

    left = aux.read.into_list(cfg['mat']['left_nodes'])
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    data = []
    
    add_data(cfg,data,['all','all'])
    add_data(cfg,data,['all','low'],edge_filter=low_pass_edge_filter,args=35)
    add_data(cfg,data,['all','mid'],edge_filter=mid_pass_edge_filter,args=(35,66))
    add_data(cfg,data,['all','high'],edge_filter=high_pass_edge_filter,args=66)
    add_data(cfg,data,['ipsilateral','all'],edge_filter=ipsilateral_pass_filter,args=[left,right])
    add_data(cfg,data,['contralateral','all'],edge_filter=contralateral_pass_filter,args=[left,right])
    
    df = pd.DataFrame(data,columns=["Comparison","Network","Edge threshold","Measure","Cell","Jaccard Distance"])
    if source_data: df.to_csv(source_data,index=False)


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
    run(params.config,source_data=SOURCE)
