"""
c4_contact_matrix.py

Overlays C4 contacts over and M4 ordered matrix

author@ Christopher Brittin
@date 30 May 2020
"""
import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
from igraph import Graph
import networkx as nx
import numpy as np
import ioaux
from collections import defaultdict 
import seaborn as sns
import scipy.cluster.hierarchy as sch
import matplotlib as mpl
from matplotlib import colors


from connectome.load import from_db
from connectome.format_graphs import *
from cluster_population_plot_figure import *
from connectome.load import reference_graphs
from networks.nxrandom import Randomize

CONFIG = 'configs/config.ini'
mpl.rcParams['ytick.labelsize'] = 6

def reorder_clusters(nodes,brainmap):
    border = ['Anterior','Lateral','Sublateral','Avoidance','Taxis']
    _nodes = []
    for cls in border:
        for n in nodes:
            if brainmap[n] == cls: _nodes.append(n)
    return _nodes
 
def build_mask(zone,size):
    [l1,l2] = zone 
    Z = np.zeros([size,size])
    
    for i in range(l1):
        z = np.ones(size-i)
        Z += np.diagflat(z,i)
        if i > 0: Z += np.diagflat(z,-1*i)
    for i in range(l1,l1+l2):
        z = 2*np.ones(size-i)
        Z += np.diagflat(z,i)
        Z += np.diagflat(z,-1*i)
    
    for i in range(l1+l2,l1+2*l2):
        z = 3*np.ones(size-i)
        Z += np.diagflat(z,i)
        Z += np.diagflat(z,-1*i)
    
    for i in range(l1+2*l2,l1+3*l2):
        z = 4*np.ones(size-i)
        Z += np.diagflat(z,i)
        Z += np.diagflat(z,-1*i)
    
    Z[np.where(Z==0)] = 5
    Z -= 1
    return Z
    
def build_matrix(A,Z,clusters):
    zorder = range(int(Z.max())+1)
    synval = int(Z.max()) + 1
    B = np.zeros(A.shape)
    zcount = np.zeros(len(zorder))
    ccount = np.zeros(len(zorder))
    for (k,z) in enumerate(zorder):
        idx = zip(*np.where(Z == z))
        for (i,j) in idx:
            if A[i,j] > 0: 
                zcount[k] += A[i,j] 
                B[i,j] = synval
                if clusters[i] == clusters[j]: 
                    if clusters[i] != 'Unclassified': ccount[k] += A[i,j]
            else:
                B[i,j] = z 

    return B,zcount,ccount

def build_source_data(A,Z,clusters,nodes):
    SD = [['pre','post','zone','synapse','intra_cluster']]
    zorder = range(int(Z.max())+1)
    for (k,z) in enumerate(zorder):
        idx = zip(*np.where(Z == z))
        for (i,j) in idx:
            is_syn = 0
            is_cluster = 0
            if A[i,j] > 0: 
                is_syn = A[i,j] 
                if clusters[i] == clusters[j]: 
                    if clusters[i] != 'Unclassified': is_cluster = A[i,j]
            SD.append([nodes[i],nodes[j],z,is_syn,is_cluster])

    return SD



def extract_multi_neigh(A,Z,clusters,nodes):
    CLUSTERS = {'Anterior':0,'Lateral':1,'Sublateral':2,'Avoidance':3,'Taxis':4,'Unclassified':5}
    pre = np.zeros([A.shape[0],6])
    post = np.zeros([A.shape[0],6])
    zorder = range(int(Z.max())+1)
    synval = int(Z.max()) + 1

    for z in [3,4]:
        idx = zip(*np.where(Z == z))
        for (i,j) in idx:
            if A[i,j] > 0: 
                ic = CLUSTERS[clusters[nodes[i]]]
                nodes[j]
                clusters[nodes[j]]
                jc = CLUSTERS[clusters[nodes[j]]]
                pre[i,jc] += A[i,j]
                post[j,ic] += A[i,j]

    mpre = sorted([nodes[i] for i in np.where(np.sum(pre,axis=1)>0)[0]])
    mpost = sorted([nodes[i] for i in np.where(np.sum(post,axis=1)>0)[0]])
    
    mneigh = list(sorted(list(set(mpre) | set(mpost))))
    return mneigh

def plot_matrix3(im,M,bundles,cfg,zones,no_cbar=False,weight='weight',vmin=0,vmax=0.3,cbar_ticks=[0,0.15,0.3]):
    nodes = sorted(M.nodes())
    idx = im.dendrogram_row.reordered_ind
    nodes = [nodes[i] for i in idx]
    bcolor = dict([(d[2],d[1]) for d in ioaux.read.into_list2(cfg['mat']['bundle_color'])])
    cclass = ioaux.read.into_dict(cfg['mat']['class'])
    #ncolor = [bcolor[bundles[n]] for n in nodes]
    brainmap = ioaux.read.into_dict(cfg['clusters']['brainmap'])
    COLORS = ['0','0.25','0.5','0.75','1','#FF69B4']#,'#959595']#'#B2B2B2']
    BOUNDS = [0,1,2,3,4,5,6]
    cmap = colors.ListedColormap(COLORS)
    norm = colors.BoundaryNorm(BOUNDS, cmap.N)
    
    bundles['CEPDL'] = 'Anterior'
    nodes = reorder_clusters(nodes,brainmap) 
    ncolor = [bcolor[bundles[n]] for n in nodes]
    ndx = dict([(n,i) for (i,n) in enumerate(nodes)])
    #A = np.zeros([len(nodes),len(nodes)])
    A = nx.to_numpy_array(M,nodelist=nodes,weight=None)
    W = nx.to_numpy_array(M,nodelist=nodes,weight='count')
    
    Z = build_mask(zones,M.number_of_nodes())
    clusters = [bundles[n] for n in nodes]
    B,zcount,ccount = build_matrix(W,Z,clusters)    
    SD = build_source_data(W,Z,clusters,nodes)    
    
    #mneigh = extract_multi_neigh(W,Z,bundles,nodes)
    #ioaux.write.from_list('source_data/source_data_figure_3abc.csv',SD)

    zsum = zcount.sum()
    csum = ccount.sum()
    ccount = np.divide(ccount,zcount) 
    zcount /= zsum
    #ccount /= csum
    print(A.sum(),W.sum())
    print(zcount,zsum)
    print(ccount,csum)
    fig,_ax = plt.subplots(2,1,figsize=(1.4,2))
    cols = ['0','0.25','0.5','0.75','1.0']
    x = [0,1,2,3,4]
    ax = _ax[0]
    ax.bar(x,zcount,color=cols,edgecolor='k')
    ax.set_ylim([0,1.])
    ax.set_ylabel('$\mathbb{C}^4$',fontsize=8)
    ax.set_yticks([0,0.25,0.5,0.75,1.0]) 
    ax.set_xticks(x)
    ax.set_xticklabels(['0','1','2','3','4'],fontsize=6)
    #ax.set_title('Fraction of synapses',fontsize=8)
    ax.text(0.4,0.8,f'$n={int(zsum)}$',transform=ax.transAxes,fontsize=6)
    
    ax = _ax[1]
    ax.bar(x,ccount,color=cols,edgecolor='k')
    ax.set_ylim([0,1.])
    ax.set_ylabel('Intra-cluster',fontsize=6)
    ax.set_yticks([0,0.25,0.5,0.75,1.0]) 
    ax.set_xticks(x)    
    ax.set_xticklabels(['0','1','2','3','4'],fontsize=6)
    ax.text(0.4,0.8,f'$n={int(csum)}$',transform=ax.transAxes,fontsize=6)



    plt.tight_layout()
    #plt.savefig('results/cluster_revision/localized_synapses_zones_zcounts.svg')
 
    #fig, ax = plt.subplots(figsize=(10,10))
    #ax.imshow(A, cmap=cmap, norm=norm)
    nodes = []
    hm = sns.clustermap(B,row_cluster=False,col_cluster=False,col_colors=ncolor,row_colors=ncolor,
            yticklabels=nodes,xticklabels=[],cbar_pos=(.06, .35, .03, .15),
            cbar_kws={"ticks":cbar_ticks},cmap=cmap
            #,figsize=(2.5,2.5))
            ,figsize=(2.5,2.5)) 
    hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize = 4)
    hm.cax.set_visible(False)
    #plt.savefig('results/cluster_revision/localized_synapses_zones.png',dpi=400)
    #plt.savefig('results/cluster_revision/localized_synapses.svg')
    


def filter_reference_graph(Ref,tid):
    H = nx.Graph()
    if Ref.is_directed(): H = nx.DiGraph()
    for (u,v,w) in Ref.edges.data(data=True):
        if w['id'] != tid: continue
        H.add_edge(u,v,weight=w['weight'],id=w['id'])
        if 'sections' in w: H[u][v]['sections'] = w['sections']
    return H


def run(params):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    clusters = ioaux.read.into_dict(params.clusters) 
    im = plot_clustermap(params.fin,cfg,clusters,
            no_cbar=params.no_cbar,no_ticklabels=params.no_ticklabels)

    A,C,E = reference_graphs(cfg)

    M = nx.read_graphml(cfg['refgraphs']['adj_cl'])
    nodes = sorted(M.nodes())
    deg = [M.degree(n) for n in nodes]
    print('mu,std',np.mean(deg),np.std(deg))
    zone = [int(0.5*np.mean(deg)) + 1, int(np.std(deg))]
    #plot_matrix(im,M,clusters,cfg,vmax=0.3)
    #plt.savefig('results/cluster_revision/m4_contacts.png',dpi=400)
    
    Ref = make_reference_graphs(C)
    Ref = make_synapse_reference_graph_dir(Ref,A[4])
    M = filter_reference_graph(Ref,4)
    print('#edges1',M.number_of_edges(),Ref.number_of_edges())
    M = collapse_lr_nodes(M,left,right)
    for m in nodes:
        if not M.has_node(m): M.add_node(m)
    print('#edges2',M.number_of_edges())
    count = 0
    for (u,v,w) in M.edges.data(data='count'): count += 2
    print('sum',count)
    #plot_matrix(im,M,clusters,cfg,weight='sections')
    #plt.show()
    rc = ['DVA', 'AIBL', 'RIBL', 'RIAL', 'AVEL', 'AVDL', 'AVAL','PVCL']
    CLUSTER = {'Anterior':0,'Lateral':1,'Sublateral':2,'Avoidance':3,'Taxis':4,'Unclassified':5}
    rr = np.zeros([8,6])
    for (i,n) in enumerate(rc):
        for m in M.neighbors(n):
            rr[i,CLUSTER[clusters[m]]] += 1
        for m in M.predecessors(n):
            rr[i,CLUSTER[clusters[m]]] += 1
    
    print('rr',rr)
 
    plot_matrix3(im,M,clusters,cfg,zone,weight='sections')
    plt.show()

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

    params = parser.parse_args()
    run(params)
