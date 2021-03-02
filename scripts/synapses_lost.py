"""
@name: synapses_lost.py 
@description:
Quantify fraction of synaptic contacts lost as function of membrane contact threshold for M4 graph.

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-05
"""
import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
import copy
import numpy as np
from collections import defaultdict

from connectome.load import from_db
from connectome.format_graphs import consensus_graph,filter_graph_edge
from connectome.format_graphs import make_reference_graphs,clean_graph
import aux

#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
FS = 7


BAD_NODES = ['PVR','PVNL','HSNL','PVR.','HSNR','SABD','VB02','PLNL','PLNR','PVNR']

def plot_synapse_lost_membrane_contact(ax,trange,syn_lost,syn_n,edge_thresh=35):
    for (k,v) in sorted(syn_lost.items()): 
        ax.plot(trange,v,linewidth=1,label='$\delta=%d$ ($n=%d$)'%(k,syn_n[k]))
    #ax.plot(trange,gap_lost,'b-',linewidth=4,label='Lost gap junctions')
    ax.legend(fontsize=6)
    ax.set_xlim([30,70])
    ax.set_xlabel('$\mathbb{M}^4$ contact area percentile',fontsize=FS)
    #ax.set_ylabel('Fraction of synapses lost',fontsize=10)
    #ax.set_title('Reference dataset $\mathbb{A}^%d$'%a_delta,fontsize=18)
    ax.grid()
    ax.set_ylim([0,1])
    ax.set_yticks(np.arange(0,1.1,0.1))
    ax.axvline(edge_thresh,linewidth=1,color='k',linestyle='--')


def create_source_data(trange,syn_lost):
    k = len(trange)
    data = np.zeros([k,5])
    data[:,0] = trange
    for i in range(1,5): data[:,i] = syn_lost[i]
    return data

def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)
    left = aux.read.into_list(cfg['mat']['left_nodes'])
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = aux.read.into_lr_dict(cfg['mat']['lrmap'])
    nodes = aux.read.into_list(cfg['mat']['nodes']) 
    remove = aux.read.into_list(cfg['mat']['remove'])
    edge_thresh = cfg.getint('params','lower_weight_threshold')
    dbs = cfg['input']['databases'].split(',')

    _G = []
    for d in dbs:
        D = from_db(d,adjacency=True,chemical=True,electrical=True,remove=remove,dataType='networkx')
        _G.append(D)

    G = []
    for _D in _G:
        D = copy.deepcopy(_D)
        D.split_left_right(left,right)
        D.map_right_graphs(lrmap)
        G.append(D)

    A,C,E = {},{},{}
    for i in range(1,5):
        A[i] = nx.Graph()
        C[i] = nx.DiGraph()
        E[i] = nx.Graph()
        consensus_graph(A[i],[G[0].Al,G[0].Ar,G[1].Al,G[1].Ar],i,nodes,weight=['weight','wnorm'])
        consensus_graph(C[i],[G[0].Cl,G[0].Cr,G[1].Cl,G[1].Cr],i,nodes)
        consensus_graph(E[i],[G[0].El,G[0].Er,G[1].El,G[1].Er],i,nodes)
    
    A4 = A[4]
    Ref = make_reference_graphs(A)
    cRef = make_reference_graphs(C)
    eRef = make_reference_graphs(E)
    
    cdelta = np.zeros(4)
    misscount = 0
    for (a,b) in cRef.edges():
        if not Ref.has_edge(a,b): 
            misscount += 1
            #print(a,b,G[0].A.has_edge(lrmap[a],lrmap[b]),G[1].A.has_edge(lrmap[a],lrmap[b]))
            continue
        aid = Ref[a][b]['id'] - 1
        cdelta[aid] += 1
    csum = cdelta.sum()
    print(cdelta)
    cdelta /= csum
    print(f'missing {misscount}')

    edelta = np.zeros(4)
    misscount = 0
    for (a,b) in eRef.edges():
        if not Ref.has_edge(a,b): 
            misscount += 1 
            continue 
        aid = Ref[a][b]['id'] - 1
        edelta[aid] += 1
    esum = edelta.sum()
    print(edelta)
    edelta /= esum
    print(f'missing {misscount}')

    fig,ax = plt.subplots(1,1,figsize=(2,2))
    ax.grid(zorder=1)
    x = np.array([1,2,3,4])
    width=0.35
    lbl = ['$\mathbb{M}^%d$'%i for i in range(1,5)] 
    ax.bar(x-width*0.5,cdelta,width=0.33,label='Syn. contacts ($n$=%d)'%int(csum),zorder=2)
    ax.bar(x+width*0.5,edelta,width=0.33,label='Gap j. contacts ($n$=%d)'%int(esum),zorder=2)
    ax.set_xlim([0.5,4.5]) 
    ax.set_xticks(x)
    ax.set_xticklabels(lbl,fontsize=FS)
    ax.set_ylabel('Fraction of synaptic contacts',fontsize=FS)
    ax.set_ylim([0,1])
    ax.legend(loc='upper left',fontsize=6)
    plt.tight_layout()
    plt.savefig('./results/consensus_graph_synapses_membrane_delta.svg')
    
    #Synaptic contacts lost in A4 graph
    trange = range(30,75,5)
    a_delta = 4
    chem_lost = defaultdict(list)
    gap_lost = defaultdict(list)
    chem_n = {}
    gap_n = {}

    for t in tqdm(trange,desc="Edge thresh"):
        G = []
        for _D in _G:
            D = copy.deepcopy(_D)   
            D.A = filter_graph_edge(D.A,pct=t)
            D.split_left_right(left,right)  
            D.map_right_graphs(lrmap)
            G.append(D)

        A = nx.Graph()
        consensus_graph(A,[G[0].Al,G[0].Ar,G[1].Al,G[1].Ar],a_delta,nodes)
        
        _n,_m = 0,0
        tmp = {1:0,2:0,3:0,4:0}
        pmp = {1:0,2:0,3:0,4:0}
        cct = {1:0,2:0,3:0,4:0}
        for i in range(1,a_delta+1):
            C = nx.DiGraph()
            E = nx.Graph()
            consensus_graph(C,[G[0].Cl,G[0].Cr,G[1].Cl,G[1].Cr],i,nodes)
            consensus_graph(E,[G[0].El,G[0].Er,G[1].El,G[1].Er],i,nodes)
            
            C = clean_graph(C,A4)
            E = clean_graph(E,A4)
           
            chem_n[i] = C.number_of_edges()
            gap_n[i] = E.number_of_edges()
            missed_syn = [(a,b,C[a][b]['weight']) for (a,b) in C.edges() if not A.has_edge(a,b)]
            #print(len(missed_syn),C.number_of_edges())
            chem_lost[i].append(float(len(missed_syn)) / C.number_of_edges())
            tmp[i] = float(len(missed_syn)) / C.number_of_edges()*100
            pmp[i] = len(missed_syn) 
            cct[i] = C.number_of_edges()
            _n += C.number_of_edges()
            _m += float(len(missed_syn))
            missed_syn = [(a,b,E[a][b]['weight']) for (a,b) in E.edges() if not A.has_edge(a,b)]
            gap_lost[i].append(float(len(missed_syn)) / E.number_of_edges())
        #print(t, _m / _n,_m,_n)#,chem_lost)
        if t == 35: print(f'M{a_delta}: synaptic contacts lost {100*_m/_n:.2f}: (C1:{tmp[1]:.2f}%, C2:{tmp[2]:.2f}%, C3:{tmp[3]:.2f}%, C4:{tmp[4]:.2f}%)')
        if t == 35: print(f'M{a_delta}: Total synaptic contacts {_n}, Total lost {_m}: (C1:{pmp[1]}/{cct[1]}, C2:{pmp[2]}/{cct[2]}, C3:{pmp[3]}/{cct[3]}, C4:{pmp[4]}/{cct[4]})')
        
        
    fig,ax = plt.subplots(1,2,figsize=(4.5,2))
    plot_synapse_lost_membrane_contact(ax[0],trange,chem_lost,chem_n) 
    ax[0].set_ylabel('Fraction of $\mathbb{C}^\delta$',fontsize=FS) 
    #ax[0].set_title('Synapses at $\mathbb{M}^4$',fontsize=FS) 
    plot_synapse_lost_membrane_contact(ax[1],trange,gap_lost,gap_n) 
    ax[1].set_ylabel('Fraction of $\mathbb{G}^\delta$',fontsize=FS) 
    #ax[1].set_title('Gap j. at $\mathbb{M}^4$',fontsize=FS)
    plt.tight_layout()
    plt.savefig('./results/consensus_graph_synapses_lost.svg')
    plt.show()
    
    header = ['membrane_contact_percentile'] + [f'delta_{i}' for i in range(1,5)]
    header = ','.join(header)
    data = create_source_data(trange,chem_lost)
    np.savetxt('source_data/source_data_extended_data_2e.csv',data,delimiter=',',header=header)
    
    data = create_source_data(trange,gap_lost)
    np.savetxt('source_data/source_data_extended_data_2f.csv',data,delimiter=',',header=header)


    
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
