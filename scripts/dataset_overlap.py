"""
@name: dataset_overlap.py
@description:
Look at overlap between datasets. 

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""

import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.labelsize'] = 7


from connectome.format_graphs import clean_graph
#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

def compute_overlap(G,H):
    g = set(G.edges()) 
    h = set(H.edges())

    u = len(g | h)
    i = len(g & h)
    return float(i) / u

def compute_intersection(G,H):
    g = set(G.edges()) 
    h = set(H.edges())
    i = len(g & h)
    return i

def compute_cons(w,Z,A4):
    w4 = clean_graph(w,A4)
    n = float(w4.number_of_edges())
    w4cons = np.zeros(5)
    for i in range(1,5):
        z = clean_graph(Z[i],A4)
        o = compute_intersection(w4,z)
        w4cons[i] = o
    
    w4cons[0] = n #np.sum(w4cons[1:])
    #w4cons /= n
    return w4cons
 
def fit_c4_cons(w):
    min_err = 10
    pmin = -1
    for p in np.arange(0,1.01,0.01):
        P = pcons(p)
        err = abs(P - w).sum()
        if err > min_err: continue
        min_err = err
        pmin = p
    return pmin,min_err


def pcons(p):
    q = 1-p
    P = np.zeros(5)
    P[0] = q**4
    P[1] = 4*p*(q**3)
    P[2] = 6*(p**2)*(q**2)
    P[3] = 4*(p**3)*(q**1)
    P[4] = p**4
    return P


def plot_counts(ax,data,ylabel=None,legend=None,width=0.15,width2=0.13):
    labels = ['$\mathbb{C}^%d$'%i for i in range(1,5)]
    ind = np.arange(4)
    ax.bar(ind-width,data[0,:],width=width2,color='r',zorder=2,label='Witvliet2020')
    ax.bar(ind,data[1,:],width=width2,color='g',zorder=2,label='White1986')
    ax.bar(ind+width,data[2,:],width=width2,zorder=2,color='b',label='Cook2019')

    #ax.set_ylim([0,0.8])
    ax.set_xlim([-0.5,3.5])
    ax.set_xticks(ind)
    if labels: ax.set_xticklabels(labels,fontsize=12)
    ax.legend(fontsize=12)
    ax.set_ylabel('Number of synapses',fontsize=12) 
    #if ylabel: ax.set_ylabel(ylabel,fontsize=6)


def plot_extended_overlap(ax,data,ylabel=None,legend=None,width=0.15,width2=0.13):
    labels = ['$\mathbb{C}^%d$'%i for i in range(5)]
    labels[0] = 'Total'
    ind = np.arange(5)
    N = np.array(data[:,0],copy=True)
    data[:,0] = data[:,1:].sum(axis=1)
    data /= N[:,None]
    ax.bar(ind-1.5*width,data[0,:],width=width2,label=r'$\mathbb{C}^1$ ($n$=%d)'%N[0])
    ax.bar(ind-0.5*width,data[1,:],width=width2,label=r'$\mathbb{C}^2$ ($n$=%d)'%N[1])
    ax.bar(ind+0.5*width,data[2,:],width=width2,label=r'$\mathbb{C}^3$ ($n$=%d)'%N[2])
    ax.bar(ind+1.5*width,data[3,:],width=width2,label=r'$\mathbb{C}^4$ ($n$=%d)'%N[3])

    ax.set_ylim([0,1.0])
    ax.set_xlim([-0.5,4.5])
    ax.set_xticks(ind)
    if labels: ax.set_xticklabels(labels,fontsize=7)
    ax.legend(fontsize=6)
    ax.set_ylabel('Number of synapses',fontsize=7) 
    #if ylabel: ax.set_ylabel(ylabel,fontsize=6)



def plot_overlap(ax,data,ylabel=None,legend=None,width=0.15,width2=0.13):
    labels = ['$\mathbb{C}^%d$'%i for i in range(1,5)]
    ind = np.arange(4)
    ax.bar(ind-0.5*width,data[0,:],width=width2,color='gold',zorder=2,label='Witvliet&White')
    ax.bar(ind+0.5*width,data[1,:],width=width2,color='darkviolet',zorder=2,label='White&Cook')

    ax.set_ylim([0,1.0])
    ax.set_xlim([-0.5,3.5])
    ax.set_xticks(ind)
    if labels: ax.set_xticklabels(labels,fontsize=12)
    ax.legend(fontsize=12)
    ax.set_ylabel('Percent overlap',fontsize=12) 
    #if ylabel: ax.set_ylabel(ylabel,fontsize=6)

def plot_fit(ax,data,fit,ylabel=None,legend=None,width=0.35,width2=0.33):
    labels = ['$\mathbb{C}^%d$'%i for i in range(5)]
    ind = np.arange(len(data))
    #ax.grid(zorder=1)
    ax.bar(ind-0.5*width,data,width=width2,color='k',zorder=2,label=legend[0])
    ax.bar(ind+0.5*width,fit,width=width2,zorder=2,color='#a20000',label=legend[1])

    #ax.set_ylim([0,0.8])
    ax.set_ylim([0,1.0])
    x = [float('0.%d'%i) for i in range(1,len(data)+1)]
    ax.set_xlim([-0.5,len(data)-0.5])
    ax.set_xticks(ind)
    if labels: ax.set_xticklabels(labels,fontsize=8)
    ax.legend(fontsize=6)
    ax.set_ylabel('Fraction of contacts',fontsize=8) 
    #if ylabel: ax.set_ylabel(ylabel,fontsize=6)



def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)
    A4 = nx.read_graphml(cfg['refgraphs']['adj']%4)
    Y,W,Z = {},{},{}
    for i in range(1,5):
        Y[i] = nx.read_graphml(cfg['refgraphs']['chem']%i)
        W[i] = nx.read_graphml(cfg['zhen']['white_chem']%i)
        Z[i] = nx.read_graphml(cfg['zhen']['zhen_chem']%i)

    ovr = np.zeros((2,4))
    count = np.zeros((3,4))

    zsum = 0
    for i in range(1,5):
        z = clean_graph(Y[i],A4)
        zsum += z.number_of_edges()
    print('sum',zsum)

  
    cookzhen = np.zeros((4,5))
    for i in range(4): cookzhen[i,:] = compute_cons(Y[i+1],Z,A4)
    missed = []
    zscored = []
    for (a,b) in Y[4].edges():
        found = 0
        for (i,z) in Z.items():
            if z.has_edge(a,b):
                found = 1
                break
        if found: 
            zscored.append((a,b))
        else:
            missed.append((a,b))
    print('Witvliet Missed')
    for (a,b) in missed:
        print(a,b,Y[4][a][b]['weight'])
    
    print('Witvliet small scored')
    smallcount = 0
    for (a,b) in zscored:
        if Y[4][a][b]['weight'] < 5: smallcount += 1
    print(f'smallcount: {smallcount}')
    
    zhencook = np.zeros((4,5))
    for i in range(4): zhencook[i,:] = compute_cons(Z[i+1],Y,A4)
    whitezhen = np.zeros((4,5))
    for i in range(4): whitezhen[i,:] = compute_cons(W[i+1],Z,A4)
    zhenwhite = np.zeros((4,5))
    for i in range(4): zhenwhite[i,:] = compute_cons(Z[i+1],W,A4)
    cookwhite = np.zeros((4,5))
    for i in range(4): cookwhite[i,:] = compute_cons(Y[i+1],W,A4)
    whitecook = np.zeros((4,5))
    for i in range(4): whitecook[i,:] = compute_cons(W[i+1],Y,A4)
    

    print(cookzhen)
    print(zhencook)
    fig,ax = plt.subplots(1,2,figsize=(4,2))
    plot_extended_overlap(ax[0],cookzhen)
    ax[0].set_title('Cook scored by Witvliet',fontsize=7)
    ax[0].set_ylabel('Cook syn. contacts',fontsize=7)
    ax[0].set_xlabel('Witvliet syn. contacts',fontsize=7)
    plot_extended_overlap(ax[1],zhencook)
    ax[1].set_title('Witvliet scored by Cook',fontsize=7)
    ax[1].set_ylabel('Witvliet syn. contacts',fontsize=7)
    ax[1].set_xlabel('Cook syn. contacts',fontsize=7)
    plt.tight_layout()
    #plt.savefig('./results/zhen_data/cook_witvliet_delta_overlap.svg')
    
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
