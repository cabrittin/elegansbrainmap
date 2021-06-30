"""
@name: cluster_perturbation_validation.py
@description:
Look at differences in weight differences to find appropriate noise level

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""
import sys
sys.path.append(r'./preprocess')
import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import networkx as nx
import numpy as np
import ioaux
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
from itertools import combinations
import multiprocessing_on_dill as mp
import time
from scipy.stats import norm

from connectome.load import from_db
from connectome.load import reference_graphs
#from connectome.format_graphs import *
from connectome.format_graphs import consensus_graph,filter_graph_edge,normalize_edge_weight
from measures import probability_dist
from cluster_population_models import *

mpl.rcParams['xtick.labelsize'] = 5 
mpl.rcParams['ytick.labelsize'] = 5 
#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

def plot_log(V,perturbed=False):    
    _V = np.log(V)
    mu = np.mean(_V[:,0])
    std = np.std(_V[:,0])
    __V = _V - mu
    _V = (_V - mu) / std

    Z = np.zeros([V.shape[0],2])
    Z[:,0] = _V[:,0]
    Z[:,1] = np.std(_V[:,1:],axis=1)
    
    
    print(mu,std)
    #mu = np.mean(Z[np.where(Z[:,1]<1)[0],1])
    mu = np.mean(Z[:,1])
    #mu = np.median(Z[:,1])
    fig,_ax = plt.subplots(1,3,figsize=(5.4,1.8))
    FS = 7
    
    ax = _ax[0]
    ax.hist(__V[:,0],bins=60,range=(-3,3),density=True)
    x = np.linspace(-5,5,100)
    y = norm(0,std)
    ax.plot(x,y.pdf(x),'r-',linewidth=1)
    ax.set_xlim([-3,3])
    ax.set_ylim([0,0.6])
    ax.set_yticks([0,0.3,0.6])
    ax.set_xlabel('log of membrane contact area',fontsize=FS)
    ax.set_ylabel('Frequency',fontsize=FS)
    ax.set_title('$\mathbb{M}^4$ contacts',fontsize=FS)
    if perturbed: ax.set_title('$\widetilde{\mathbb{M}}^4$ contacts',fontsize=FS)
    ax.text(0.5,0.8,'STD$=%1.2f$'%std,transform=ax.transAxes,fontsize=6)

    ax = _ax[1]  
    ax.scatter(Z[:,0],Z[:,1],s=1)
    ax.set_ylim([0,2])
    ax.set_yticks([0,1,2])
    #ax.axhline(1,color='r',linestyle='dashed')
    ax.set_xlabel('log of mean',fontsize=FS)
    ax.set_ylabel('Standard deviation',fontsize=FS) 
    ax.set_title('Variability across 4 datasets',fontsize=FS)

    ax = _ax[2]
    ax.hist(Z[:,1],bins=40,range=(0,2))
    ax.set_ylim([0,150])
    ax.set_yticks([0,50,100,150])
    ax.set_ylabel("# of contacts sites",fontsize=FS)
    ax.set_xlabel("Standard deviation",fontsize=FS)
    ax.axvline(mu,linestyle='dashed',color='r')
    ax.text(0.4,0.8,"Mean$=%1.2f$"%mu,transform=ax.transAxes,fontsize=6)
    #plt.savefig('results/cluster_revision/noise_distribution.svg')
    plt.tight_layout()
    #plt.show()

def plot_linear(V):

    Z = np.zeros([V.shape[0],2])
    Z[:,0] = V[:,0]
    Z[:,1] = np.std(V[:,1:],axis=1)
    
    #idx = np.where(Z[:,0] < 40000)[0]
    #Z = Z[idx,:]
    #p = np.polyfit(Z[:,0],Z[:,1],1)
    
    #print(p)
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    x = np.linspace(0,Z.max(),10) 
    ax.scatter(Z[:,0],Z[:,1]/Z[:,0],s=2)
    #ax.plot(x,x,'r-')
    #ax.plot(x,x*p[0] + p[1],'g-')
    #mx = Z.max()
    ax.set_xlim([0,Z[:,0].max()])
    #ax.set_ylim([0,Z[:,1].max()])
    ax.set_xlabel('Mean contact area')
    ax.set_ylabel('Std/mean contact area')
    plt.tight_layout() 
    #plt.savefig('results/cluster_revision/linear_quotient_noise.svg')
    
    #mu = np.mean(Z[np.where(Z[:,1]<1.2)])
    #mu = np.mean(Z[:,1])
    #fig,_ax = plt.subplots(1,3,figsize=(7.5,2))
    #FS = 10
    
def plot_noise(_ax,V):
    vmax = np.max(V[:,1:])
    ax = _ax[0] 
    for i in range(2,5):
        ax.scatter(V[:,1]/vmax,V[:,i]/vmax,s=2,color='b')
    z = np.linspace(0,1,11)
    ax.plot(z,z,'r-')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_xlabel('Equivalent contact',fontsize=10)
    ax.set_ylabel('L4 left',fontsize=10)
    ax.set_title('Membrane contact areas')
 
    ax = _ax[1]
    for i in range(2,5):
        y = np.log(np.true_divide(V[:,1],V[:,i]))
        x = V[:,1] / np.max(V[:,1])
        ax.scatter(x,y,s=2,color='b')
    ax.axhline(0,color='r',linewidth=1)
    ax.set_xlim([0,1])
    ax.set_ylim([-2,2])
    ax.set_xlabel('L4 left',fontsize=10)
    ax.set_ylabel('log(Equivalent/L4 left)',fontsize=10)
    ax.set_title('Membrane contact areas')
    plt.tight_layout() 
 
def run(_cfg,fout=None,source_data=None):
    wrn = """NOTE: This script pulls directly from the sql databases.
    If these databases are not setup then this script will crash"""
    
    print(wrn)
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)
    
    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = ioaux.read.into_lr_dict(cfg['mat']['lrmap'])
    nodes = ioaux.read.into_list(cfg['mat']['nodes']) 
    remove = ioaux.read.into_list(cfg['mat']['remove'])
    edge_thresh = cfg.getint('params','lower_weight_threshold')
    dbs = cfg['input']['databases'].split(',')

    G = []
    for d in dbs:
        D = from_db(d,adjacency=True,chemical=True,electrical=True,remove=remove,dataType='networkx')
        D.A = filter_graph_edge(D.A,pct=edge_thresh)
        D.split_left_right(left,right)  
        D.map_right_graphs(lrmap)
        G.append(D)
    
    H = [G[0].Al,G[0].Ar,G[1].Al,G[1].Ar]
    for g in H: normalize_edge_weight(g)
    M = nx.Graph() 
    consensus_graph(M,H,4,nodes,weight=['weight','wnorm'])
    sd = []
    V = np.zeros([M.number_of_edges(),5])
    for (i,(u,v,w)) in enumerate(M.edges.data('wnorm')):
        V[i,0] = w
        tmp = [0,0,0,0]
        for (j,h) in enumerate(H): 
            V[i,j+1] = h[u][v]['wnorm']
            tmp[j] = h[u][v]['wnorm']
        sd.append([u,v,w] + tmp)
    
    ioaux.write.from_list('source_data/ed_fig5_emperical.csv',sd)
    
    lscale = get_log_scale(cfg)
    sig = 0.23
    _gsizes,Mp,Hp = perturb_data(cfg,['N2U','JSH'],lscale,sig,spatial_domain=0)
    sd = []
    Vp = np.zeros([Mp.number_of_edges(),5])
    for (i,(u,v,w)) in enumerate(Mp.edges.data('wnorm')):
        Vp[i,0] = w
        tmp = [0,0,0,0]
        for (j,h) in enumerate(Hp): 
            Vp[i,j+1] = h[u][v]['wnorm']
            tmp[j] = h[u][v]['wnorm']
        sd.append([u,v,w] + tmp)
    
    ioaux.write.from_list('source_data/ed_fig5_pert.csv',sd)
    
    plot_log(V)
    #plt.savefig('results/cluster_revision/noise_distribution.svg')
    plot_log(Vp,perturbed=True)
    #plt.savefig('results/cluster_revision/noise_distribution_perturbed.svg')
   
    
    #fig,_ax = plt.subplots(1,2,figsize=(5,2))
    #plot_noise(_ax,V)
    
    #fig,_ax = plt.subplots(1,2,figsize=(5,2))
    #plot_noise(_ax,Vp)

    fig,_ax = plt.subplots(1,2,figsize=(3.6,1.8))
    x = np.zeros([M.number_of_edges(),2])
    for (i,(u,v,w)) in enumerate(M.edges.data('wnorm')):
        x[i,0] = w
        if Mp.has_edge(u,v): x[i,1] = Mp[u][v]['wnorm']
    
    ax = _ax[0]
    ax.scatter(x[:,0]/np.max(x),x[:,1]/np.max(x),s=1)
    z = np.linspace(0,1,11)
    ax.plot(z,z,'r-')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_xlabel('Empirical',fontsize=7)
    ax.set_ylabel('Perturbed',fontsize=7)
    ax.set_title('Membrane contact areas',fontsize=7)
     
    
    ax = _ax[1]
    x[:,1] = np.log(np.true_divide(x[:,1],x[:,0]))
    ax.scatter(x[:,0]/np.max(x[:,0]),x[:,1],s=1)
    ax.axhline(0,color='r',linewidth=1)
    ax.set_xlim([0,1])
    ax.set_ylim([-2,2])
    ax.set_xlabel('Empirical',fontsize=7)
    ax.set_ylabel('log(Perturbed/Empirical)',fontsize=7)
    ax.set_title('Membrane contact areas',fontsize=7)
    plt.tight_layout() 
    #plt.savefig('results/cluster_revision/empirical_vs_perturbed.svg')
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
