"""
cluster_population_plot_figure.py

Monte Carlo community analysis

author@ Christopher Brittin
@date 30 May 2019
"""
import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
from igraph import Graph
import networkx as nx
import numpy as np
import aux
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
from connectome.format_graphs import consensus_graph,filter_graph_edge


CONFIG = 'configs/config.ini'
mpl.rcParams['ytick.labelsize'] = 6

def reorder_clusters(nodes,brainmap):
    border = ['Anterior','Lateral','Sublateral','Avoidance','Taxis']
    _nodes = []
    for cls in border:
        for n in nodes:
            if brainmap[n] == cls: _nodes.append(n)
    return _nodes

def plot_clustermap(fin,cfg,bundles,no_ticklabels=False,no_cbar=False,reorder=False):
    fout = fin.replace('.npz','.png')
    norder = np.load(fin,allow_pickle=True)['norder']
    c = np.load(fin,allow_pickle=True)['C'] #[1000 x n] array of perturbed clusters
    gsizes = np.load(fin,allow_pickle=True)['gsizes']

    norder = dict([(i,n) for (n,i) in dict(norder.tolist()).items()])
    nodes = [norder[i] for i in range(len(norder))] 

    n = c.shape[1]
    z = np.zeros([n,n])
    
    for (i,j) in combinations(range(n),2):
        s = (c[:,i] == c[:,j]).sum()
        z[i,j] = s
        z[j,i] = s
    m = nx.read_graphml(cfg['refgraphs']['adj_cl'])
    nodes = sorted(m.nodes())
    cclass = aux.read.into_dict(cfg['mat']['class'])
    #bundles = aux.read.into_dict(cfg['clusters']['noise'])
    bcolor = dict([(d[2],d[1]) for d in aux.read.into_list2(cfg['mat']['bundle_color'])])
    border = aux.read.into_list(cfg['mat']['bundle_order'])
    ngroup = defaultdict(list)
    for (k,v) in bundles.items(): ngroup[v].append(k)
    for k in ngroup: ngroup[k] = sorted(ngroup[k])
    norder = []
    for k in border: norder += ngroup[k]

    z /= c.shape[0]
    np.fill_diagonal(z,1)
    y = sch.linkage(z, method='ward')#,optimal_ordering = True)
    d = sch.dendrogram(y, orientation='right',no_plot=True)
    ncolor = [bcolor[bundles[n]] for n in nodes]
    _nodes = [cclass[n] for n in nodes]
    if no_ticklabels: nodes = [] 
    im = sns.clustermap(z,row_linkage=y,col_linkage=y,dendrogram_ratio=(.2, .2),
            cbar_pos=(.12, .82, .03, .15),xticklabels=[],yticklabels=_nodes,
            #row_colors=ncolor,col_colors=ncolor,figsize=(2.5,2.5))
            row_colors=ncolor,col_colors=ncolor,figsize=(10,10))
    im.ax_row_dendrogram.set_visible(False)
    if no_cbar:
        im.cax.set_visible(False)
    else: 
        im.ax_cbar.tick_params(labelsize=10)
        im.ax_cbar.set_label('cluster frequency')
    plt.savefig(fout,dpi=400)
    im.reordered_ind = im.dendrogram_row.reordered_ind

    if reorder:
        idx = im.dendrogram_row.reordered_ind
        nodes = [nodes[i] for i in idx]
        nodes = reorder_clusters(nodes,bundles) 
        #ndx = dict([(n,i) for (i,n) in enumerate(nodes)])
        ndx = dict([(n,i) for (i,n) in enumerate(nodes)])
        ncolor = [bcolor[bundles[n]] for n in nodes]
        _nodes = [cclass[n] for n in nodes]
        _z = np.zeros(np.shape(z))
        NODES = sorted(m.nodes())
        for (i,u) in enumerate(NODES):
            for (j,v) in enumerate(NODES):
                _z[ndx[u],ndx[v]] = z[i,j]
        _y = sch.linkage(_z, method='single')
        _im = sns.clustermap(_z,row_cluster=False,col_cluster=False,
                #row_linkage=_y,col_linkage=_y,dendrogram_ratio=(.2, .2),
                cbar_pos=(.05, .9, .20, .03),xticklabels=[],yticklabels=[],#_nodes,
                row_colors=ncolor,col_colors=ncolor,figsize=(10,10),
                cbar_kws={"orientation": "horizontal","ticks":[0,0.2,0.4,0.6,0.8,1.0]})
        _im.ax_cbar.tick_params(labelsize=10)
        axx = _im.ax_col_dendrogram.axes
        axx.clear()
        link = im.dendrogram_col.linkage
        _link = link[:,:]
        #np.savetxt('scratch/temp.csv',_link,delimiter=',')
        _link[89,:] = np.array([181,179,18.6392725848366,46])
        _link[90,:] = np.array([182,180,23.9234768351246,70])
        #link[[4, 2]] = link[[2, 4]] 
        #print(link[[4,2]])
        #print(link[[2,4]])
        with plt.rc_context({'lines.linewidth': 1.0}):
            sch.dendrogram(_link,ax=axx,orientation='top',link_color_func=lambda x: 'k')
        axx.set_yticklabels(['']*len(axx.get_yticklabels()))
        axx.tick_params(color='w')
        ndx = dict([(n,i) for (i,n) in enumerate(NODES)])
        im.reordered_ind = [ndx[n] for n in nodes]
        fout = fin.replace('.npz','_reorder.png')
        plt.savefig(fout,dpi=400)

    return im

def plot_compare(cfg,im):
    COLORS = ['#9301E7', '#E7A401', '#5E7FF1','#FC0000','#1FFF00','#9b9b9b','#199b07']
    CLUSTER = {"Anterior":0,"Lateral":1,"Sublateral":2,"Avoidance":3,
            "Taxis":4,"Unclassified":5,"Taxis2":6}
    _CLUSTER = dict([(v,k) for (k,v) in CLUSTER.items()])
    BOUNDS = [0,1,2,3,4,5,6,7]
    cmap = colors.ListedColormap(COLORS)
    norm = colors.BoundaryNorm(BOUNDS, cmap.N)
    
    L4 = aux.read.into_dict('data/clusters/clusters_s23_JSH_t35.csv')
    M4 = aux.read.into_dict('data/clusters/clusters_s23_m4_t35.csv')
    Adult = aux.read.into_dict('data/clusters/clusters_s23_N2U_t35.csv')
    
    #L4 = aux.read.into_dict('data/clusters/clusters_s23_JSH_t35_opt.csv')
    #M4 = aux.read.into_dict('data/clusters/clusters_s23_m4_t35_opt.csv')
    #Adult = aux.read.into_dict('data/clusters/clusters_s23_N2U_t35_opt.csv')


    #idx = im.dendrogram_row.reordered_ind
    idx = im.reordered_ind
    m = nx.read_graphml(cfg['refgraphs']['adj_cl'])
    nodes = sorted(m.nodes())
    cclass = aux.read.into_dict(cfg['mat']['class'])
    
    bcolor = dict([(d[2],d[1]) for d in aux.read.into_list2(cfg['mat']['bundle_color'])])
    nodes = [nodes[i] for i in idx]
    ncolor = [COLORS[CLUSTER[Adult[n]]] for n in nodes]

    N = len(nodes)
    X = np.zeros([N,4])
    data = []
    for (i,n) in enumerate(nodes):
        X[i,0] = CLUSTER[M4[n]]
        X[i,1] = CLUSTER[L4[n]]
        X[i,2] = CLUSTER[Adult[n]]
        X[i,3] = 5
        if np.min(X[i,:3]) == np.max(X[i,:3]): X[i,3] = np.max(X[i,:3])
        data.append([i,cclass[n],M4[n],L4[n],Adult[n],_CLUSTER[X[i,3]]])
    #aux.write.from_list('results/cluster_revision/fig1_row_order.csv',data)
    #aux.write.from_list('results/cluster_revision/edfig5i_row_order.csv',data)
    
    fig, ax = plt.subplots(figsize=(3.8,9))
    ax.imshow(X, cmap=cmap, norm=norm)
    ax.set_aspect('auto')
    ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
    ax.set_xticks([0.5,1.5,2.5,3.5])
    #ax.set_xticklabels(['Perturbed contacts','$\mathbb{M}^4$','L4','Adult'],fontsize=10,ha='left')
    ax.set_xticklabels([])
    ax.set_yticks([])
    plt.text(0.17,0.07,'$\widetilde{\mathbb{M}}^4$', fontsize=18, transform=plt.gcf().transFigure)
    plt.text(0.37,0.07,'$\widetilde{L4}$', fontsize=18, transform=plt.gcf().transFigure)
    plt.text(0.50,0.07,'$\widetilde{Adult}$', fontsize=18, transform=plt.gcf().transFigure)
    plt.text(0.73,0.07,'Final', fontsize=18, transform=plt.gcf().transFigure)
    #plt.text(0.77,0.07,'Final', fontsize=18, transform=plt.gcf().transFigure)
    #plt.savefig('results/cluster_revision/cluster_compare.svg')
 

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('fin',
                        action = 'store',
                        help = 'Input file: perturbation npz file')
    
    parser.add_argument('clusters',
                        action = 'store',
                        help = 'Cluster file')
  
    
    parser.add_argument('-c','--config',
                    dest = 'config',
                    action = 'store',
                    default = CONFIG,
                    required = False,
                    help = 'Config file')
    
    parser.add_argument('--no_cbar',
                    dest = 'no_cbar',
                    action = 'store_true',
                    default = False,
                    required = False,
                    help = 'Remove colorbar')
    
    parser.add_argument('--no_ticklabels',
                    dest = 'no_ticklabels',
                    action = 'store_true',
                    default = False,
                    required = False,
                    help = 'Remove cell ticklabels')
    
    parser.add_argument('--reorder_nodes',
                    dest = 'reorder',
                    action = 'store_true',
                    default = False,
                    required = False,
                    help = 'Reorder nodes')



    params = parser.parse_args()
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    clusters = aux.read.into_dict(params.clusters) 
    im = plot_clustermap(params.fin,cfg,clusters,
            no_cbar=params.no_cbar,no_ticklabels=params.no_ticklabels,reorder=params.reorder)

    plot_compare(cfg,im) 
    
    plt.show()
    
