"""
@name: breakdown_contact_sites.py
@description:
Looks at the distribution of contact sites. 

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""

import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from connectome.load import reference_graphs
from connectome.format_graphs import *
from measures import probability_dist

mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'
FS = 7

def plot_consensus_degree_dist(ax,adj,title=None,xlim=[-2,2],ylim=[0,1],label='M'):
    col = ['tab:purple','tab:blue','tab:orange','tab:green','tab:red']
    px,x = probability_dist(adj[:,1])        
    source_data = np.zeros([len(x),6])
    ax.plot(x,px,linestyle='--',linewidth=1,color='k',
                label=r'All')
    source_data[:,0] = x
    source_data[:,5] = px
    imin = int(adj[:,0].min())
    imax = int(adj[:,0].max()) 
    for i in range(1,imax+1):
        _adj = adj[adj[:,0] == i]
        px,_ = probability_dist(_adj[:,1])
        source_data[:,i] = px
        ax.plot(x,px,linewidth=1,color=col[i],label=r'$\mathbb{%s}^{%d}$'%(label,i))

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend(loc='upper right',fontsize=6)
    ax.set_ylabel('Survival distribution',fontsize=FS)
    ax.set_xlabel('Standardized $\mathbb{M}^4$ contact area',fontsize=FS)
    if title: ax.set_title(title,fontsize=7)
    return source_data
    
def plot_edge_counts(ax,ecount,width=0.2):
    x = np.arange(4)
    ax.grid(zorder=0)
    esum = ecount.sum(axis=1)
    ax.bar(x-width,ecount[0,:],width=width,label=r'$\mathbb{M}^{\delta}$ ($n$=%d)'%esum[0],zorder=3)
    ax.bar(x,ecount[1,:],width=width,label=r'$\mathbb{C}^{\delta}$ ($n$=%d)'%esum[1],zorder=3)
    ax.bar(x+width,ecount[2,:],width=width,label=r'$\mathbb{G}^{\delta}$ ($n$=%d)'%esum[2],zorder=3)
    ax.set_ylabel('# of contacts',fontsize=FS)
    ax.legend(fontsize=6,loc='upper left')
    ax.set_xticklabels(['$\delta=$%d' %d for d in range(1,5)],fontsize=7)
    ax.set_xticks(x)
 
def filter_reference_graph(Ref,tid):
    H = nx.Graph()
    if Ref.is_directed(): H = nx.DiGraph()
    for (u,v,w) in Ref.edges.data(data=True):
        if w['id'] != tid: continue
        H.add_edge(u,v,weight=w['weight'],id=w['id'])
        if 'sections' in w: H[u][v]['sections'] = w['sections']
    return H


def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)
    
    A,C,E = reference_graphs(cfg)
    
    ecounts = np.zeros((3,4))
    for i in range(1,5):
        j = i - 1
        ecounts[0,j] = A[i].number_of_edges()
        c = clean_graph(C[i],A[4])
        ecounts[1,j] = c.number_of_edges()
        e = clean_graph(E[i],A[4])
        ecounts[2,j] = e.number_of_edges()
    
    print(ecounts)

    fig,ax = plt.subplots(1,1,figsize=(2.1,2))
    plot_edge_counts(ax,ecounts)
    plt.tight_layout()
    plt.savefig('./results/contact_dist.svg')
    header = 'memb_contact_area,delta_1,delta_2,delta_3,delta_4,all'
    Ref = make_reference_graphs(A)
    mu = np.mean([w for (a,b,w) in Ref.edges.data('weight')])
    print(mu,mu*450*1e-6)
    standardize_edge_weigth(Ref)
    adj = graph_to_array(Ref,attr=['id','weight']) 
    px,x = probability_dist(adj[:,1])
    fig,ax = plt.subplots(1,1,figsize=(2,2))
    source_data = plot_consensus_degree_dist(ax,adj)
    ax.set_xlabel('Standardized $\mathbb{M}^\delta$ contact area',fontsize=FS)
    plt.tight_layout()
    plt.savefig('./results/contact_adj_membrane_survival.svg') 
    np.savetxt('./source_data/source_data_extended_data_2h.csv',source_data,delimiter=',',header = header)

    Ref = make_reference_graphs(C)
    Ref = make_synapse_reference_graph_dir(Ref,A[4])
    #print('_M edges',_M.number_of_edges())
    mu = np.mean([w for (a,b,w) in Ref.edges.data('weight')])
    print(mu,mu*450*1e-6)
    standardize_edge_weigth(Ref)
    adj = graph_to_array(Ref,attr=['id','weight']) 
    px,x = probability_dist(adj[:,1])
    fig,ax = plt.subplots(1,1,figsize=(2,2))
    source_data = plot_consensus_degree_dist(ax,adj,label='C')
    plt.tight_layout()
    plt.savefig('./results/contact_chem_membrane_survival.svg') 
    np.savetxt('./source_data/source_data_extended_data_2i.csv',source_data,delimiter=',',header = header)
   
    Ref = make_reference_graphs(E)
    Ref = make_synapse_reference_graph_dir(Ref,A[4])
    mu = np.mean([w for (a,b,w) in Ref.edges.data('weight')])
    print(mu,mu*450*1e-6)
    standardize_edge_weigth(Ref)
    adj = graph_to_array(Ref,attr=['id','weight']) 
    px,x = probability_dist(adj[:,1])
    fig,ax = plt.subplots(1,1,figsize=(2,2))
    source_data = plot_consensus_degree_dist(ax,adj,label='G')
    plt.tight_layout()
    plt.savefig('./results/contact_gap_membrane_survival.svg') 
    np.savetxt('./source_data/source_data_extended_data_2j.csv',source_data,delimiter=',',header = header)
    plt.show()

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

    run(params.config)
