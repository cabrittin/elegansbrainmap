"""
mc_community_synapse.py

Monte Carlo community analysis

author@ Christopher Brittin
@date 30 May 2019
"""
import sys
sys.path.append(r'./analysis')
import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
from igraph import Graph
import networkx as nx
import numpy as np
import ioaux
from collections import defaultdict 
from sklearn.metrics.cluster import normalized_mutual_info_score
import matplotlib.pyplot as plt
from tqdm import tqdm
from itertools import combinations
import seaborn as sns
import scipy.cluster.hierarchy as sch
import matplotlib as mpl
from matplotlib import colors


from connectome.load import from_db
from connectome.format_graphs import *
from connectome.load import reference_graphs
from cluster_population_plot_figure import *

#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'
mpl.rcParams['ytick.labelsize'] = 6
    
def plot_matrix(im,M,bundles,cfg,no_cbar=False,weight='weight',vmin=0,vmax=0.3,
        fig_title=None,cbar_ticks=[0,0.15,0.3]):
    nodes = sorted(M.nodes())
    idx = im.dendrogram_row.reordered_ind
    nodes = [nodes[i] for i in idx]
    bcolor = dict([(d[2],d[1]) for d in ioaux.read.into_list2(cfg['mat']['bundle_color'])])
    ncolor = [bcolor[bundles[n]] for n in nodes]
    brainmap = ioaux.read.into_dict(cfg['clusters']['brainmap'])
    nodes = reorder_clusters(nodes,brainmap)
    ncolor = [bcolor[bundles[n]] for n in nodes]
    A = nx.to_numpy_array(M,nodelist=nodes,weight=weight)
    A = A/A.sum(axis=1)[:,None]
    A[np.isnan(A)] = 0

    nodes = []
    hm = sns.clustermap(A,row_cluster=False,col_cluster=False,col_colors=ncolor,row_colors=ncolor,
            yticklabels=nodes,xticklabels=[],cbar_pos=(.06, .35, .03, .15),
            cbar_kws={"ticks":cbar_ticks},
            vmin=vmin,vmax=vmax,figsize=(2.5,2.5))
    if no_cbar: hm.cax.set_visible(False)
    if fig_title: hm.fig.canvas.set_window_title(fig_title)

def filter_reference_graph(Ref,tid):
    H = nx.Graph()
    if Ref.is_directed(): H = nx.DiGraph()
    for (u,v,w) in Ref.edges.data(data=True):
        if w['id'] != tid: continue
        H.add_edge(u,v,weight=w['weight'],id=w['id'])
        if 'sections' in w: H[u][v]['sections'] = w['sections']
    return H

def reorder_clusters(nodes,brainmap):
    border = ['Anterior','Lateral','Sublateral','Avoidance','Taxis']
    _nodes = []
    for cls in border:
        for n in nodes:
            if brainmap[n] == cls: _nodes.append(n)
    return _nodes
 
def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)

    perturbations = 'data/perturbations/mc_cluster_rand_sig23_m4_t35.npz' 
    #clusters = 'data/clusters/clusters_s23_m4_t35.csv'
    clusters = 'data/clusters/final_clusters.csv'
    no_ticklabels = False
    reorder_nodes = True
    no_cbar = False
 
    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    clusters = ioaux.read.into_dict(clusters) 
    im = plot_clustermap(perturbations,cfg,clusters,
            no_cbar=no_cbar,no_ticklabels=no_ticklabels)

    A,C,E = reference_graphs(cfg)

    M = nx.read_graphml(cfg['refgraphs']['adj_cl'])
    nodes = sorted(M.nodes())
    plot_matrix(im,M,clusters,cfg,vmax=0.3,fig_title='ED 7a M4')
    #plt.savefig('results/cluster_revision/m4_contacts.png',dpi=400)
    
    M = nx.read_graphml(cfg['refgraphs']['animal_adj_cl']%'JSH')
    plot_matrix(im,M,clusters,cfg,no_cbar=True,vmax=0.3,fig_title='ED 7a L4')
    #plt.savefig('results/cluster_revision/jsh_contacts.png',dpi=400)
    
    M = nx.read_graphml(cfg['refgraphs']['animal_adj_cl']%'N2U')
    plot_matrix(im,M,clusters,cfg,no_cbar=True,vmax=0.3,fig_title='ED 7a Adult')
    #plt.savefig('results/cluster_revision/n2u_contacts.png',dpi=400)
    
    Ref = make_reference_graphs(C)
    Ref = make_synapse_reference_graph_dir(Ref,A[4])
    
    for i in range(1,5):
        M = filter_reference_graph(Ref,i)
        M = collapse_lr_nodes(M,left,right)
        for m in nodes:
            if not M.has_node(m): M.add_node(m)
        _cbar = True
        if i == 4: _cbar = False
        plot_matrix(im,M,clusters,cfg,no_cbar=_cbar,weight='sections',fig_title='ED 7b C%d'%i)
        #plt.savefig('results/cluster_revision/m4_c%d_contacts.png'%i,dpi=400)
    
    Ref = make_reference_graphs(E)
    Ref = make_synapse_reference_graph_dir(Ref,A[4])
    
    for i in range(1,5):
        M = filter_reference_graph(Ref,i)
        M = collapse_lr_nodes(M,left,right)
        for m in nodes:
            if not M.has_node(m): M.add_node(m)
        _cbar = True
        if i == 4: _cbar = False
        plot_matrix(im,M,clusters,cfg,no_cbar=_cbar,weight='sections',fig_title='ED 7c G%d'%i)
        #plt.savefig('results/cluster_revision/m4_g%d_contacts.png'%i,dpi=400)
        
    plt.show()

if __name__=="__main__":
    run(CONFIG) 
