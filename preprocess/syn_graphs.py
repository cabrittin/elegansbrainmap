"""
@name: syn_graphs.py
@description:
    format synaptic graphs 

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-04
"""
import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
import networkx as nx
from collections import defaultdict

from expand_list_to_bilateral import expand_list
from connectome.format_graphs import *
import aux


#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

def clean_graph(G,A):
    edges = [(a,b) for (a,b) in G.edges()]
    for (a,b) in edges:
        if not A.has_edge(a,b):
            G.remove_edge(a,b)

def graph_v1(G,bundle):
    H = nx.DiGraph()
    for (u,v) in G.edges():
        if G.node[u]['bundle'] == bundle: 
            intra = int(G.node[u]['bundle'] == G.node[v]['bundle'])
            if not intra: continue
            H.add_edge(u,v,weight=float(G[u][v]['weight']),
                        intra=int(G[u][v]['intra']),
                        bundle=G[u][v]['bundle'])
    return H


def score_nodes(G,thresh=0.05):
    #bunin = defaultdict(float)
    #bunout = defaultdict(float)
    bundeg = defaultdict(float)
    for n in G.nodes():
        bundeg[G.nodes[n]['merge']] += G.in_degree(n,weight='weight')
        bundeg[G.nodes[n]['merge']] += G.out_degree(n,weight='weight')
        bundeg['deg'] += (G.in_degree(n,weight='weight') + G.out_degree(n,weight='weight'))
     
    for n in G.nodes():
        cin = G.in_degree(n,weight='weight')
        cout = G.out_degree(n,weight='weight')
        G.nodes[n]['deg'] = min(0.3,(cin + cout) / bundeg[G.nodes[n]['merge']])
        G.nodes[n]['fan'] = 'in'
        if cin == 0:
            G.nodes[n]['fan'] = 'neutral'
            if  cout > 0: G.nodes[n]['fan'] = 'out'
            continue
        if cout > 0:
            fancoef = np.log(cin / cout)


lbl = ['bundle','merge','axon','class','cook']
bndl = ['Anterior','Taxis','Avoidance','Lateral','Sublateral']
#DEG = 3
rc = ['AVAL','AVBL','AVDL','AVEL','PVCL','DVA']
rc2 = ['AIAL','AIBL','AIZL','AVAL','AVBL','AVEL','RIAL','RICL','RIML','RIPL','RMDVL','SMDVL']
rm = ['PVR','HSNL','HSNR','PLNL','PLNR','PVNL','PVNR','PVR.']

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')

    parser.add_argument('-d','--deg',
                        dest = 'deg',
                        action = 'store',
                        type = int,
                        default = 4,
                        required = False,
                        help = 'Synapse degree of reproducibility')
    
    parser.add_argument('-m','--m_delta',
                        dest = 'mdelta',
                        action = 'store',
                        type = int,
                        default = 4,
                        required = False,
                        help = 'Membrane contact degree of reproducibility')
    
    parser.add_argument('--cell',
                        dest = 'cell',
                        action = 'store',
                        default = None,
                        required = False,
                        help = 'Cell name')
    
    parser.add_argument('--l3',
                        dest = 'l3',
                        action = 'store_true',
                        default = False,
                        required = False,
                        help = 'Assess layer 3 directionality')
    
    parser.add_argument('--rc',
                        dest = 'rc',
                        action = 'store_true',
                        default = False,
                        required = False,
                        help = 'Assess rich club neurons')
    
    parser.add_argument('--rc2',
                        dest = 'rc2',
                        action = 'store_true',
                        default = False,
                        required = False,
                        help = 'Assess rich club neurons 2')
 
    params = parser.parse_args()
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
   
    #bundles = aux.read.into_list2(cfg['bundles']['bundles'])
    bundles = aux.read.into_list2(cfg['clusters']['final'])
    merge = aux.read.into_list2(cfg['clusters']['brainmap'])
    axclass = aux.read.into_list2(cfg['bundles']['class'])
    left = aux.read.into_list(cfg['mat']['left_nodes'])
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    rlmap = dict(zip(right,left))
    lrdict = aux.read.into_lr_dict(cfg['mat']['lrmap'])
    nclass = aux.read.into_dict(cfg['mat']['class'])
    cook = aux.read.into_dict(cfg['mat']['cook'])
    resnet = aux.read.into_dict(cfg['mat']['resnet'])

    bundles = expand_list(bundles,lrdict)
    merge = expand_list(merge,lrdict)
    axclass = expand_list(axclass,lrdict)
    
    master = {}
    for (i,[n,b]) in enumerate(bundles):
        tmp = [b,merge[i][1],axclass[i][1],nclass[n],cook[n]]
        master[n] = dict(zip(lbl,tmp))

    params.deg = min(params.mdelta,params.deg) 
    A = nx.read_graphml(cfg['refgraphs']['adj']%params.mdelta)
    C = nx.read_graphml(cfg['refgraphs']['chem']%params.deg)
    #nx.write_weighted_edgelist(C,'./mat/c4_edgelist.csv',delimiter=',')
    A.remove_nodes_from(rm)
    C.remove_nodes_from(rm)
    if params.rc:
        C = nx.DiGraph()
        for i in range(3,5):
            C = nx.compose(C,nx.read_graphml(cfg['refgraphs']['chem']%i))

    clean_graph(C,A)
    #nx.add_node_attributes(C,master)
    
    #Read resnet class
    redges = aux.read.into_list2(cfg['resnet']['resnet_edges'])
    R = nx.DiGraph()
    for (u,v,c) in redges:
        if u in rlmap: u = rlmap[u]
        if v in rlmap: v = rlmap[v] 
        R.add_edge(u,v,rtype=c)

    H = nx.DiGraph()
    for (u,v) in C.edges():
        w = C[u][v]['weight']
        if u in rlmap: u = rlmap[u]
        if v in rlmap: v = rlmap[v]
        if not H.has_edge(u,v): H.add_edge(u,v,weight=0)
        H[u][v]['weight'] += w
    nx.set_node_attributes(H,master)
    try:
        for (u,v) in H.edges():
            H[u][v]['intra'] = int(H.nodes[u]['merge'] == H.nodes[v]['merge'])
            H[u][v]['bundle'] = H.nodes[v]['bundle']
            H[u][v]['axon'] = H.nodes[v]['axon']
            layer3 = int(int(resnet[u]) == 3 and int(resnet[v]) == 3)
            H[u][v]['in_resnet'] = int(layer3 or H[u][v]['intra'])
            H[u][v]['brainmap'] = int(H[u][v]['intra'] and (resnet[u]==resnet[v]))
            direction = 'resnet'
            if not H[u][v]['in_resnet']:
                direction = 'v_ff'
                if int(resnet[v]) < int(resnet[u]): direction = 'v_fb'
            H[u][v]['direction'] = direction
            H[u][v]['rtype'] = 'None'
            if R.has_edge(u,v): 
                H[u][v]['rtype'] = R[u][v]['rtype']
            else:
                print(u,v)

    except Exception as e:
        print(e)
        print('Not using M4 data, therefore skipping edge assignments')
   
    H.remove_edges_from(nx.selfloop_edges(H))
    score_nodes(H)

    data = aux.read.into_list2(cfg['motifs']['classified'])
    M = nx.DiGraph()
    mcells = []
    sc = defaultdict(int)
    tc = defaultdict(int)
    ic = defaultdict(int)
    for d in data:
        [a,b,c] = d[:3]
        if a in right: a = lrdict[a]
        if b in right: b = lrdict[b]
        if c in right: c = lrdict[c]
        sc[a] += 1
        tc[c] += 1
        ic[b] += 1
        M.add_edge(a,b)
        M.add_edge(b,c)
        M.add_edge(a,c)
        mcells += [a,b,c]
    """ 
    print(f'mcells: {len(set(mcells))}')
    print(np.median(list(sc.values())),np.median(list(ic.values())),np.median(list(tc.values())))
    print({k: v for k, v in sorted(sc.items(), key=lambda item: item[1])})
    print({k: v for k, v in sorted(ic.items(), key=lambda item: item[1])})
    print({k: v for k, v in sorted(tc.items(), key=lambda item: item[1])})
    print(', '.join(sorted([nclass[k] for (k,v) in sc.items() if v > 2])))
    print(', '.join(sorted([nclass[k] for (k,v) in ic.items() if v > 2])))
    print(', '.join(sorted([nclass[k] for (k,v) in tc.items() if v > 2])))
    """
    for (a,b) in H.edges():
        H[a][b]['in_motif'] = 0
        if M.has_edge(a,b): H[a][b]['in_motif'] = 1
    
    in_mot = 0.
    num,den = 0.,0.
    for (a,b,w) in H.edges.data(data=True): 
        if H[a][b]['in_resnet'] > 0:
            den += 1
            if H[a][b]['in_motif'] > 0: num += 1
    #print('% in motif',in_mot,H.number_of_edges(),in_mot / H.number_of_edges(),M.number_of_edges())
    print('% in motif',num,den,num/den)
    fout = cfg['syngraphs']['chem']%(params.mdelta,params.deg)
    if params.cell:
        fout = fout.replace('.graphml',f'_{params.cell}.graphml') 
        h = [params.cell] + [n for n in H.neighbors(params.cell)] + [n for n in H.predecessors(params.cell)]
        H = H.subgraph(h)
    elif params.l3:
        l3 = [n for n in H.nodes() if int(resnet[n]) == 3]
        fout = cfg['syngraphs']['chem']%params.deg
        fout = fout.replace('.graphml',f'_l3.graphml')
        H = H.subgraph(l3)
        L = nx.DiGraph()
        for (u,v) in H.edges():
            if not H[u][v]['in_motif']: continue
            bu = H.node[u]['merge']
            bv = H.node[v]['merge']
            L.add_edge(bu,bv,weight=H[u][v]['weight'])
        for n in L.nodes(): 
            L.node[n]['bundle'] = n
            L.node[n]['class'] = n
        H = L 
    
    elif params.rc:
        fout = cfg['syngraphs']['chem']%params.deg
        fout = fout.replace('.graphml',f'_rc.graphml')
        H = H.subgraph(rc) 
    elif params.rc2:
        fout = cfg['syngraphs']['chem']%(params.mdelta,params.deg)
        fout = fout.replace('.graphml',f'_rc2.graphml')
        H = H.subgraph(rc2) 
   
    #u = 'OLLL'
    #for v in H.neighbors(u):
    #    print(u,v,H[u][v]['in_motif'])
    print(f"Write to: {fout}")
    nx.write_graphml(H,fout)
         
