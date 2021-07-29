"""
@name: format_graphs.py
@description:

Module for formating and manipulating graphs

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-03
"""

import networkx as nx
import numpy as np

def consensus_graph_deprecated(G,H,deg,nodes):
    """
    Make degree consensus consensus graphs. 
    
    Looks for edges found in only n=deg number of H graphs

    Input:
    ------
    G : Graph to store consensus edges
    H : List of networkx graphs
    deg: int, number of graphs required to have edge
    nodes: list of nodes to build consensus graphs

    """
    for n in nodes:
        neigh = set([])
        for h in H:
            if not h.has_node(n): continue
            neigh = neigh.union(set(list(h.neighbors(n))))
        for m in neigh:
            _deg = 0
            _w = 0.
            for h in H:
                if h.has_edge(n,m):
                    _deg += 1
                    _w += h[n][m]['weight']
            if _deg == deg:
                G.add_edge(n,m,weight=_w/_deg)
        if not G.is_directed(): continue
        for h in H:
            if not h.has_node(n): continue
            neigh = neigh.union(set(list(h.predecessors(n))))
        for m in neigh:
            _deg = 0
            _w = 0.
            for h in H:
                if h.has_edge(m,n):
                    _deg += 1
                    _w += h[m][n]['weight']
            if _deg == deg:
                G.add_edge(m,n,weight=_w/_deg)

def consensus_graph(G,H,_deg,nodes,weight=['weight']):
    """
    Make degree consensus consensus graphs. 
    
    Looks for edges found in only n=deg number of H graphs

    Input:
    ------
    G : Graph to store consensus edges
    H : List of networkx graphs
    deg: int, number of graphs required to have edge
    nodes: list of nodes to build consensus graphs

    """
    n = len(nodes)
    k = len(H)
    ndict = dict([(i,j) for (i,j) in enumerate(nodes)])
    Z = np.zeros([k,n,n])
    for w in weight: 
        for (i,h) in enumerate(H): 
            Z[i,:,:] = nx.to_numpy_array(h,nodelist=nodes,weight=w)
        deg = np.count_nonzero(Z,axis=0)
        S = Z.sum(axis=0)
        nz = np.nonzero(deg)
        S[nz] = np.divide(S[nz],deg[nz])
        for (i,j) in zip(*np.where(deg == _deg)):
            if not G.has_edge(nodes[i],nodes[j]): G.add_edge(nodes[i],nodes[j])
            G[nodes[i]][nodes[j]][w] = S[i,j]
        
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
        um = rlmap[u]
        vm = rlmap[v]
        if not H.has_edge(um,vm): H.add_edge(um,vm,count=0)
        H[um][vm]['count'] += 1
        for (k,w) in G[u][v].items():
            try:
                H[um][vm][k] += w
            except:
                H[um][vm][k] = w
    return H


def low_pass_edge_filter(G,pct,attr='weight'):
    """
    Remove edges with attr above therehold
    
    Inputs:
    -------
    G : igraph, directed or undirected
    pct: percentile threshold, float, between 0 and 1
    attr: str, edge attribute (deftaul 'weight')
    """
    weights = [e[attr] for e in G.es()]
    thresh = np.percentile(weights,pct)
    G.delete_edges([e for e in G.es() if e[attr] > thresh])

def mid_pass_edge_filter(G,pct):
    """
    Remove edges wiith attr between low and high threshold
    
    Inputs:
    -------
    G : igraph, directed or undirected
    pct: percentile threshold, list, [low threshold, high threshold]
    attr: str, edge attribute (deftaul 'weight')
    """
    weights = [e['weight'] for e in G.es()]
    thresh1 = np.percentile(weights,pct[0])
    thresh2 = np.percentile(weights,pct[1])
    edges = [e for e in G.es() if (e['weight'] >= thresh1 and e['weight'] <= thresh2) ]
    G.delete_edges(edges)

def high_pass_edge_filter(G,pct):
    """
    Remove edges with attr above therehold
    
    Inputs:
    -------
    G : igraph, directed or undirected
    pct: percentile threshold, float, between 0 and 1
    attr: str, edge attribute (deftaul 'weight')
    """
    weights = [e['weight'] for e in G.es()]
    thresh = np.percentile(weights,pct)
    G.delete_edges([e for e in G.es() if e['weight'] < thresh])

def low_pass_edge_filter_nx(G,pct,attr='weight'):
    """
    Remove edges with attr above therehold
    
    Inputs:
    -------
    G : networkx, directed or undirected
    pct: percentile threshold, float, between 0 and 1
    attr: str, edge attribute (deftaul 'weight')
    """
    weights = [w for (a,b,w) in G.edges.data(attr)]
    thresh = np.percentile(weights,pct)
    G.remove_edges_from([(a,b) for (a,b,w) in G.edges.data(attr) if w >= thresh])

def mid_pass_edge_filter_nx(G,pct,attr='weight'):
    """
    Remove edges wiith attr between low and high threshold
    
    Inputs:
    -------
    G : networkx, directed or undirected
    pct: percentile threshold, list, [low threshold, high threshold]
    attr: str, edge attribute (deftaul 'weight')
    """
    weights = [w for (a,b,w) in G.edges.data(attr)]
    thresh1 = np.percentile(weights,pct[0])
    thresh2 = np.percentile(weights,pct[1])
    edges = [(a,b) for (a,b,w) in G.edges.data(attr) if (w < thresh1 or w >= thresh2) ]
    print(G.number_of_edges(),len(edges))
    G.remove_edges_from(edges)

def high_pass_edge_filter_nx(G,pct,attr='weight'):
    """
    Remove edges with attr above therehold
    
    Inputs:
    -------
    G : networkx, directed or undirected
    pct: percentile threshold, float, between 0 and 1
    attr: str, edge attribute (deftaul 'weight')
    """
    weights = [w for (a,b,w) in G.edges.data(attr)]
    thresh = np.percentile(weights,pct)
    G.remove_edges_from([(a,b) for (a,b,w) in G.edges.data(attr) if w < thresh])

def filter_graph_edge(A,pct=50,thresh_high=True):
    """
    Remove edges from graphs that do satify edge threshold

    Input:
    ------
    A : input networkx graph
    pct: int, edge percentile threshold
    thresh_high: bool (default=True), if true pct is the upper threshold,
        if false, pct is the lower threshold

    Return:
    -------
    H : modified graph
    """
    H = nx.Graph()
    nodes = sorted(A.nodes())
    weights = [A[u][v]['weight'] for (u,v) in A.edges()]
    thresh = np.percentile(weights,pct)
    #thresh = 1109 
    for (u,v) in A.edges():
        if thresh_high: 
            c = thresh > A[u][v]['weight']
        else:
            c = thresh <= A[u][v]['weight']
        if c: continue
        H.add_edge(u,v)
        for (a,b) in A[u][v].items(): H[u][v][a] = b
    return H

def normalize_edge_weight(A):
    """
    Normalize edge weights by the total weight in the graph
    
    Adds edge attribute 'wnorm' to the graph

    Input:
    ------
    A : input networkx graph
    
    """
    tot = np.sum([A[u][v]['weight'] for (u,v) in A.edges()])
    for (u,v) in A.edges(): A[u][v]['wnorm'] = A[u][v]['weight'] / tot
 
def standardize_edge_weigth(G,attr='weight'):
    """
    Standardizes the edge weight of teh graph

    Log-normalizes the edge weights and stardardizes to set mean to x = 0

    Input:
    -----
    G: input networkx graph

    """
    #for e in G.edges.data('weight'): print(e)
    mu = np.mean([w for (a,b,w) in G.edges.data('weight')])
    print(f"Mean weight: {mu}  square microns: {mu*450*1e-6}")
    weight = [np.log(w) for (a,b,w) in G.edges(data=attr)]
    mu = np.mean(weight)
    std = np.std(weight)
    for (a,b,w) in G.edges(data=attr):
        G[a][b][attr] = (np.log(w) - mu) / std
 
def graph_to_array(G,attr=['weight']):
    """
    Convert graph to edge array

    Input:
    ------
    G: input graphs
    attr: list of attributes

    Return:
    -------
    arr : [n,3] numpy array where n is the number of edges, row format
    
    """
    arr = np.zeros((G.number_of_edges(),len(attr)))
    for (i,(a,b)) in enumerate(G.edges()):
        for (j,atr) in enumerate(attr): 
            arr[i,j] = G[a][b][atr]
    return arr

def make_reference_graphs(G,remove=[]):
    """
    Inputs:
    -------
    G: dictionary of networkx graphs. Note all
        grapsh are assumed to be either directed or undirected
    
    remove: list of edges to exclude from the reference graph
    """
    for (i,g) in G.items(): is_dir = g.is_directed()
    if is_dir:
        H = nx.DiGraph()
    else:
        H = nx.Graph()
    
    for (i,g) in G.items():
        for (a,b) in g.edges():
            if (a,b) in remove or (b,a) in remove: continue
            H.add_edge(a,b,weight=g[a][b]['weight'],id=i)

    return H

def make_weak_adjacency_graph(A,A4,remove=[],_cutoff=2):
    """
    Maker reference graph
    Inputs:
    -------
    A: An already formatted reference graph
    """
    d = nx.networkx.all_pairs_shortest_path_length(A4,cutoff=_cutoff)
    for (a,l) in d:
        for (b,p) in l.items():
            if (a,b) in remove or (b,a) in remove: continue
            if not A.has_edge(a,b): A.add_edge(a,b,weight=-3,id=0)
    return A


def make_synapse_reference_graph_udir(S,A4):
    """
    Maker reference graph
    Inputs:
    -------
    S: An already formatted reference graph
    A4: A4 adjacency graph use to prune the space

    """
    H = nx.Graph()
    for (a,b) in A4.edges():
        if S.has_edge(a,b):
            H.add_edge(a,b,weight=A4[a][b]['weight'],id=S[a][b]['id'],sections=S[a][b]['weight'])
        else:
            H.add_edge(a,b,weight=A4[a][b]['weight'],id=0,sections=0)
    return H


def make_synapse_reference_graph_dir(S,A4):
    """
    Maker reference graph
    Inputs:
    -------
    S: An already formatted reference graph
    A4: A4 adjacency graph use to prune the space

    """
    H = nx.DiGraph()
    for (a,b) in A4.edges():
        if S.has_edge(a,b):
            H.add_edge(a,b,weight=A4[a][b]['weight'],id=S[a][b]['id'],sections=S[a][b]['weight'])
        else:
            H.add_edge(a,b,weight=A4[a][b]['weight'],id=0,sections=0)
        if S.has_edge(b,a):
            H.add_edge(b,a,weight=A4[a][b]['weight'],id=S[b][a]['id'],sections=S[b][a]['weight'])
        else:
            H.add_edge(b,a,weight=A4[a][b]['weight'],id=0,sections=0)

    return H

def clean_graph(G,Ref):
    """
    Returns graph that has edges in both G and Ref

    Inputs:
    -------
    G : networkx, directed graph
    Ref: networkx, graph
    
    Return:
    -------
    H : networkx graph, same type as G with node and edge attributes 
    """
    H = G.copy()
    H.remove_nodes_from([n for n in G if n not in Ref])
    H.remove_edges_from([(a,b) for (a,b) in G.edges() if not Ref.has_edge(a,b)])
    return H


def clean_ref_graph(Ref,A):
    """
    Removed edges in Ref graph not in A

    Parameters:
    -----------
    Ref: networkx Graph, reference graph
    A: networkx Graph, graph for screening edges

    Return:
    -------
    Cleaned copy of Ref graph
    """
    H = nx.Graph()
    if Ref.is_directed(): H = nx.DiGraph()
    for (a,b) in Ref.edges():
        if A.has_edge(a,b):
            H.add_edge(a,b,weight=Ref[a][b]['weight'],id=Ref[a][b]['id'])
    return H


