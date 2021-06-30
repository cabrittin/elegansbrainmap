"""
mc_community.py

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
import ioaux
import matplotlib.pyplot as plt
from tqdm import tqdm
from itertools import combinations
import multiprocessing_on_dill as mp
import time

#from connectome.load import from_db
import connectome.load 
from connectome.load import reference_graphs
#from connectome.format_graphs import *
from connectome.format_graphs import consensus_graph,filter_graph_edge,normalize_edge_weight
from measures import probability_dist



#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'


def multilevel_community(fin):
    H = Graph.Read_GraphML(fin)
    vc = H.community_multilevel(weights='wnorm')
    H.vs['cluster'] = vc.membership
    d = dict([(v['id'],int(v['cluster'])) for v in H.vs])
    return d
    
def collapse_lr_nodes(G,left,right):
    rlmap = dict(zip(right,left))
    if G.is_directed():
        H = nx.DiGraph()
    else:
        H = nx.Graph()

    nodes = set(G.nodes())
    lr = set(left + right)
    dif = nodes - lr
    
    for d in dif: rlmap[d] = d
    for l in left: rlmap[l] = l
    
    for (u,v) in G.edges:
        for (k,w) in G[u][v].items():
            um = rlmap[u]
            vm = rlmap[v]
            if not H.has_edge(um,vm): H.add_edge(um,vm)
            try:
                H[um][vm][k] += w
            except:
                H[um][vm][k] = w
    return H

def apply_mask(Mask,G):
    rm_edges = [(a,b) for (a,b) in G.edges() if not Mask.has_edge(a,b)]
    G.remove_edges_from(rm_edges)
    return G

def run_mc(idx,niters,cfg,dbs,sig,spatial_domain=0,delta=-1,mask=None):
    np.random.seed(idx*np.random.randint(10000)) 
    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    gfile = f'data/gtmp{idx}.graphml'
    M = nx.read_graphml(cfg['refgraphs']['adj_cl'])
    fout = f'data/perturbations/mc_cluster_rand_{idx}.npz'
    if mask: F = nx.read_graphml(mask)
    #First: get the logscaling factors
    #Already determined that this mehtod yields reasonable values
    lscale =  get_log_scale(cfg,dbs)
    
    norder = dict([(n,i) for (i,n) in enumerate(sorted(M.nodes()))])
    C = np.zeros([niters,M.number_of_nodes()])
    
    gsizes = []
    for i in tqdm(range(niters),desc=f"{idx} iters: "):
        _gsizes,M,Hp = perturb_data(cfg,dbs,lscale,sig,spatial_domain=spatial_domain,delta=delta)
        #print('edges',M.number_of_edges())
        gsizes.append(_gsizes)
        if mask: M = apply_mask(F,M)
        M = collapse_lr_nodes(M,left,right)     
        nx.write_graphml(M,gfile)
        cls = multilevel_community(gfile)
        for (node,comm) in cls.items(): C[i,norder[node]] = comm
  
    np.savez(fout,norder=norder,idx=idx,C=C,gsizes=np.array(gsizes))
 
def perturb_data(cfg,dbs,lscale,sig,spatial_domain=0,delta=-1):
    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = ioaux.read.into_lr_dict(cfg['mat']['lrmap'])
    nodes = ioaux.read.into_list(cfg['mat']['nodes'])
    remove = ioaux.read.into_list(cfg['mat']['remove'])
    edge_thresh = cfg.getint('params','lower_weight_threshold')
    #dbs = cfg['input']['databases'].split(',')
    DEG = len(dbs)*2
    G = []
    for d in dbs:
        loader_method = getattr(connectome.load, cfg[d]['load'])
        if cfg[d]['load'] == 'from_graphml': d = cfg[d]
        D = loader_method(d,adjacency=True,chemical=False,electrical=False,
                remove=remove,dataType='networkx',spatial_domain=spatial_domain)
        for (u,v) in D.A.edges(): 
            D.A[u][v]['weight'] *= np.exp(np.random.normal(scale=sig)*lscale)
        D.A = filter_graph_edge(D.A,pct=edge_thresh)
        D.A.remove_nodes_from(['PLNL','PVR','SABD'])
        D.split_left_right(left,right)  
        D.map_right_graphs(lrmap)
        G.append(D)

    M = [nx.Graph() for i in range(DEG)]
    gsizes = []
    H = [G[0].Al,G[0].Ar]
    if DEG == 4: H = [G[0].Al,G[0].Ar,G[1].Al,G[1].Ar]

    for (i,m) in enumerate(M):
        deg = i + 1
        for g in G: 
            normalize_edge_weight(g.Al)
            normalize_edge_weight(g.Ar)
            consensus_graph(m,H,deg,nodes,weight=['weight','wnorm'])
        m.remove_nodes_from(remove)
        gsizes.append(m.number_of_edges()) 
    
    if delta < 0: delta = DEG - 1
    #print('delta',delta)
    return gsizes,M[delta],H

def get_log_scale(cfg,dbs,lower_log_thresh=4):
    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    remove = ioaux.read.into_list(cfg['mat']['remove'])
    edge_thresh = cfg.getint('params','lower_weight_threshold')
    dbs = cfg['input']['databases'].split(',')

    G = []
    for d in dbs:
        loader_method = getattr(connectome.load, cfg[d]['load'])
        if cfg[d]['load'] == 'from_graphml': d = cfg[d]
        D = loader_method(d,adjacency=True,chemical=False,electrical=False,
                remove=remove,dataType='networkx')
        D.A = filter_graph_edge(D.A,pct=edge_thresh)
        D.split_left_right(left,right)  
        G.append(D)
 
    H = [G[0].Al,G[0].Ar,G[1].Al,G[1].Ar]
    mu,std = [],[]
    for i in range(4):
        w = np.array([w for (u,v,w) in H[i].edges.data('weight')])
        w = np.log(w)
        idx = np.where(w > lower_log_thresh)
        _mu,_std = np.mean(w[idx]),np.std(w[idx])
        mu.append(_mu)
        std.append(_std)
    return np.mean(std)
 
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-c','--config',
            dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('-d','--dbs',
                        dest = 'db',
                        action = 'store',
                        default = 'JSH,N2U',
                        required = False,
                        help = 'Databases, e.g. db1,db2,db3...')
    
    parser.add_argument('-n','--nproc',
                        dest = 'nproc',
                        action = 'store',
                        default = 10,
                        type=int,
                        required = False,
                        help = 'Number of parallel processes')
    
    parser.add_argument('-i','--iter',
                        dest = 'iter',
                        action = 'store',
                        default = 100,
                        type=int,
                        required = False,
                        help = 'Number of iterations per process')
    
    parser.add_argument('-s','--sigma',
                        dest = 'sig',
                        action = 'store',
                        default = '0.23',
                        required = False,
                        help = 'Noise values, e.g. sig1,sig2,sig3....')
    
    parser.add_argument('--spatial_domain',
                        dest = 'spatial_domain',
                        action = 'store',
                        default = 0,
                        type=int,
                        required = False,
                        help = 'Spatial domain: 0->All, 1->NR, 2->VG')
    
    parser.add_argument('--mask',
                        dest = 'mask',
                        action = 'store',
                        default = None,
                        required = False,
                        help = ('Path to graphml file. If supplied, '
                                'clustering is only done for edges in the mask graphml')
                        )

    parser.add_argument('-o','--fout',
                        dest = 'fout',
                        action = 'store',
                        default = None,
                        required = False,
                        help = 'Output npz file')
 

    params = parser.parse_args()
    dbs = params.db.split(',')
    tag = 'm4'
    if len(dbs) == 1: tag = dbs[0]
    SIG = map(float,params.sig.split(','))
    niters = params.iter
    nproc = params.nproc 
    delta = -1
    #tag = 'm%d'%(delta+1)

    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    gfile = 'data/gtmp.graphml'
    print(f"Clustering datasets: {params.db}")
    print(f"Spatial domain: {params.spatial_domain}")
    fin = 'data/perturbations/mc_cluster_rand_%d.npz'
    fout = 'data/perturbations/mc_cluster_rand_sig%d_%s.npz'
    if params.spatial_domain > 0:
        fout = 'data/perturbations/mc_cluster_rand_sig%d_%s_s' + str(params.spatial_domain) + '.npz'

    deg = len(dbs)*2
    for sig in SIG:
        print(f'Sigma = {sig}')
        procs = []
        for idx in range(nproc):
            proc = mp.Process(target=run_mc, args=(idx,niters,cfg,dbs,sig,params.spatial_domain,delta,params.mask))
            procs.append(proc)
            proc.start()
        
        for proc in procs: proc.join() 
        
        
        C = np.zeros([nproc*niters,93])
        g = np.zeros([nproc*niters,deg])
        for idx in range(nproc):
            _start = idx * niters
            _end = (idx + 1) * niters
            D = np.load(fin%idx, allow_pickle=True)
            C[_start:_end,:] = D['C']
            g[_start:_end,:] = D['gsizes']
            norder = D['norder']
        
        np.savez(fout%(int(sig*100),tag),C=C,gsizes=g,norder=norder)
        print("Writting to %s"%(fout%(int(sig*100),tag)))
        for idx in range(nproc): os.remove(fin%idx)
