"""
@name: measures.py
@description:
Module to store functions to make measurements

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-05
"""

import numpy as np

def bundle_membrane_contacts(Z,A,bundles):
    """
    Breakdown membrane contact across bundles
    
    Parameters:
    -----------
    Z : (n,n) numpy array to hold the data
    A : networkx graph, directed or undirected
    bundles : dictionary for bundles format (k,v) = (cell,integer id)
    """
    for (u,v) in A.edges():
        uid = bundles[u]
        vid = bundles[v]
        w = A[u][v]['weight']
        Z[uid,vid] += w
        if not A.is_directed(): Z[vid,uid] += w

def cell_bundle_membrane_contacts(cell,A,bname,bundles):
    n = len(set(bname.values()))
    z = np.zeros(n)
    neigh = []
    for m in A.neighbors(cell):
        bid = bundles[m]
        z[bid] += A[cell][m]['weight']
        neigh.append([bid,bname[m],m,A[cell][m]['weight']])
    return z,neigh

def display_bundle_contact(neigh,s):
    rfrac = 1. / len(neigh)
    for r in sorted(neigh):
        r[3] /= s
        above = 37
        if r[3] >= rfrac: above = 31
        print("\033[1;%d;40m %s -- %s -- %2.3f"%(above,r[1],r[2],r[3]))
 
def probability_dist(data,N=None,xrange_=(-3,3),dx=0.1):
    if not N: N = len(data)
    N = float(N)
    x = np.arange(xrange_[0],xrange_[1],dx)
    px = [len(data[data>_x])/N for _x in x]
    return np.array(px),x
 
