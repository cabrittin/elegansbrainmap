"""
@name: localization_prob.py
@description:
Compute spatial reproduciblity of contacts

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""

import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from collections import defaultdict
from mpl_toolkits.axes_grid1 import make_axes_locatable

import aux
from connectome.format_graphs import make_reference_graphs
from localization import *

CONFIG = os.environ['CONFIG']


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('cell',
                        action='store',
                        help='Cell name')

    parser.add_argument('-c','--config',
            dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('-n','--nbins',
                        dest = 'nbins',
                        action = 'store',
                        type = int, 
                        default = 51,
                        required = False,
                        help = 'Number of bins')
    
    params = parser.parse_args()
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    left = aux.read.into_list(cfg['mat']['left_nodes'])
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = aux.read.into_lr_dict(cfg['mat']['lrmap'])
    remove = aux.read.into_list(cfg['mat']['remove'])
    cclass = aux.read.into_dict(cfg['mat']['class'])
    
    A4 = nx.read_graphml(cfg['refgraphs']['adj']%4)
    data = aux.read.into_list2(cfg['adj_align']['fout'])

    S = dict([(i,nx.read_graphml(cfg['refgraphs']['chem']%i)) for i in range(1,5)])
    Ref = make_reference_graphs(S)
    syn = aux.read.into_list2(cfg['adj_align']['cout']) 
    syn4 = format_syn_data(syn,Ref,rfilter=4) 
    syn1 = format_syn_data(syn,Ref,rfilter=3)#,high_pass=0) 
    
    #S = dict([(i,nx.read_graphml(cfg['refgraphs']['gap']%i)) for i in range(1,5)])
    #Ref = make_reference_graphs(S)
    #syn = aux.read.into_list2(cfg['adj_align']['gout']) 
    #gap4 = format_syn_data(syn,Ref,thresh=4) 
    #gap1 = format_syn_data(syn,Ref,thresh=1,high_pass=0) 

    edict = {}
    for (i,(a,b)) in enumerate(A4.edges()):
        edict[i] = (a,b)
        A4[a][b]['idx'] = i
        
    D = build_loc_data(params.cell,data,A4,nbins=params.nbins) 
    S4 = build_syn_data(params.cell,syn4,A4,nbins=params.nbins)
    S1 = build_syn_data(params.cell,syn1,A4,nbins=params.nbins)
    #G4 = build_syn_data(params.cell,gap4,A4)
    #G1 = build_syn_data(params.cell,gap1,A4)
    neigh = sorted([v for v in A4.neighbors(params.cell)])
    f = prob_response(D)
    fs4 = prob_response(S4)
    fs1 = prob_response(S1)
    #fg4 = prob_response(G4)
    #fg1 = prob_response(G1)
    
    #F = row_response(D)
    #O = breakdown_overlap(D)
    #S = format_sequences(D,neigh)
    #for s in S: print(s)
    #lcs = S[0]
    #for i in range(1,4):
    #    L = build_lcs(lcs,S[i])
    #    lcs = extract_lcs(lcs,S[i],L)
    #print(lcs,L[-1][-1])
    #print(np.mean(F,axis=0),np.std(F,axis=0))
    #label = ['L4 left','L4 right','Adult left','Adult right']
    
    """ 
    fig,ax = plt.subplots(4,2,figsize=(10,10)) 
    plot_raster(ax[:,0],D,label=label) 
    ax[3,0].set_xlabel('Effective $z$',fontsize=10) 
    ax[0,1].axis('off') plot_contact_rate(ax[1,1],f) 
    ax[1,1].set_title(cclass[params.cell]) 
    plot_contact_dist(ax[2,1],F) 
    #plot_order(ax[0,1],lcs,D[0].shape)   
    plot_order2(ax[3,1],D)   
    plt.tight_layout() 
    plt.savefig('results/%s_localization_lcs.png'%params.cell) 
    #fig,ax = plt.subplots(1,1,figsize=(6,6)) 
    #plot_order2(ax,D) 
    #plt.savefig(f'results/{params.cell}_raster.svg') 
    """ 
    fig,ax = plt.subplots(1,2,figsize=(8,2.5))
    Z = plot_order2(ax[0],D)
    ax[0].set_title(f'Localization of {cclass[params.cell]} membrane contacts')
    plot_contact_rate2(ax[1],f)
    plt.tight_layout()
    plt.savefig(f'results/adj_loc_syn/{params.cell}_adj_nbins{params.nbins}.svg')

    fig,ax = plt.subplots(1,2,figsize=(8,2.5))
    plot_order2(ax[0],S4,mask=Z)
    ax[0].set_title('Localization of %s $\mathbb{C}^4$ synapses'%cclass[params.cell])
    plot_contact_rate2(ax[1],fs4)
    plt.tight_layout()
    plt.savefig(f'results/adj_loc_syn/{params.cell}_c4_nbins{params.nbins}.svg')
    
    #fig,ax = plt.subplots(1,2,figsize=(8,2.5))
    #plot_order2(ax[0],S1,mask=Z)
    #ax[0].set_title('Localization of %s $\mathbb{C}^1$ synapses'%params.cell)
    #plot_contact_rate2(ax[1],fs1)
    #plt.tight_layout()
    #plt.savefig(f'results/adj_loc_syn/{params.cell}_c3_nbins{params.nbins}.svg')
    
    """
    fig,ax = plt.subplots(1,2,figsize=(8,2.5))
    plot_order2(ax[0],G4)
    ax[0].set_title('Localization of %s $\mathbb{G}^4$ synapses'%params.cell)
    plot_contact_rate2(ax[1],fg4)
    plt.tight_layout()
    
    fig,ax = plt.subplots(1,2,figsize=(8,2.5))
    plot_order2(ax[0],G1)
    ax[0].set_title('Localization of %s $\mathbb{G}^1$ synapses'%params.cell)
    plot_contact_rate2(ax[1],fg1)
    plt.tight_layout()
    """ 

    plt.show()

