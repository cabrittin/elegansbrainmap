"""
hierarchy.py

Module for computing hierarchy scores.

created: Christopher Brittin
date: 14 November 2018

"""
import numpy as np
import networkx as nx

def hierarchy(G,root,thresh,max_levels=15):
    root = set(root)
    H = {}
    for n in G.nodes():
        print(n)
        if n in root:
            H[n] = 0
        else:
            H[n] = -1
    for l in range(max_levels):
        H = _hierarchy(l,H,G,thresh)
    
    return H
            
def _hierarchy(l0,H,G,thresh):
    root = set([n for n in H.keys() if H[n] <= l0 and H[n] >= 0])
    print(root)
    l = l0 + 1
    for n in H.keys():
        if H[n] > -1: continue
        neigh = set(G.neighbors(n))
        nodes = list(neigh and root)
        if not nodes: continue
        s,sl = 0,0
        for m in neigh: 
            s += G[n][m]['weight']
            if m in nodes: sl += G[n][m]['weight']
        ratio = float(sl)/s
        if ratio >= thresh:
            #print n,s,sl,float(sl)/s
            H[n] = l    
    return H


def varshney_modified(H,G):
    """
    Modified version of the Varshney hierarchy 
    
    https://doi.org/10.1371/journal.pcbi.1001066

)    Parameters:
    ----------
    H : dcict
     Hierarchy scores (from hierarchy)
    G : Networkx
     Graph object

    Returns:
    --------
    nodes : list 
      Ordered list of node names. Order corresponds to xyz
    xyz : numpy array [N,3] where N in the number of nodes
      (x,y,z) coordinates
    
    """

    N = G.number_of_nodes()
    nodes = sorted(G.nodes())
    nk = dict([(nodes[i],i) for i in range(N)])
    z = np.zeros(N)
    
    for i in range(N):
        z[i] = H[nodes[i]]

    W = nx.to_numpy_matrix(G,nodelist=nodes,weight='weight')
    
    d = np.sum(W,axis=1)
    D = np.multiply(np.eye(N),d)
    
    L = D - W

    D2 = np.sqrt(D)
    D2 = np.linalg.pinv(D2)

    Q = np.dot(D2,np.dot(L,D2))

    w,v = np.linalg.eig(Q)

    idx = np.argsort(w)
    w = w[idx]
    v = v[idx]

    xyz = np.zeros([N,4])
    xyz[:,0] = v[1]
    xyz[:,1] = v[2]
    xyz[:,2] = z

    for i in sorted(list(set(H.values()))):
        layer = [n for n in H.keys() if H[n] == i]
        idx = np.array([nk[n] for n in layer])
        x = xyz[idx,1]
        jdx = np.argsort(x)
        idx = idx[jdx]
        M = len(idx)
        if M > 1:
            xlin = np.linspace(-1,1,M)
            for i in range(M):
                xyz[idx[i],3] = xlin[i]
        else:
            xyz[idx[0],3] = 0

    return nodes,xyz
