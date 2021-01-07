"""
@name: classification_measures 

@description
Functions used for precision and specificity measures

    Data structure key:
    True positive = 1
    True Negative = 2
    False positive = 3
    False negative = 4

@author: Christopher Brittin
@email: 'cabrittin' + <at> + 'gmail' + '.' + 'com'
"""

import numpy as np
import networkx as nx

def probability_dist(data,N=None,xrange_=(-3,3),dx=0.1):
    if not N: N = len(data)
    N = float(N)
    x = np.arange(xrange_[0],xrange_[1],dx)
    px = [len(data[data>_x])/N for _x in x]
    return np.array(px),x

def is_true_positive(a,b,G,Ref):
    flag = 0
    c1 = Ref[a][b]['id'] > 2
    c2 = G.has_edge(a,b)
    if c1 and c2: flag = 1
    return flag

def is_true_negative(a,b,G,Ref):
    flag = 0
    c1 = Ref[a][b]['id'] < 3
    c2 = not G.has_edge(a,b)
    if c1 and c2: flag = 1
    return flag

def is_false_positive(a,b,G,Ref):
    flag = 0
    c1 = Ref[a][b]['id'] < 3
    c2 = G.has_edge(a,b)
    if c1 and c2: flag = 1
    return flag

def is_false_negative(a,b,G,Ref):
    flag = 0
    c1 = Ref[a][b]['id'] > 2
    c2 = not G.has_edge(a,b)
    if c1 and c2: flag = 1
    return flag

def format_measures(G,Ref):
    n = Ref.number_of_edges()
    M = np.zeros((n,2))

    for (i,(a,b)) in enumerate(Ref.edges()):
        w = Ref[a][b]['weight']
        if is_true_positive(a,b,G,Ref):
            M[i,:] = [1,w]
        elif is_true_negative(a,b,G,Ref):
            M[i,:] = [2,w]
        elif is_false_positive(a,b,G,Ref):
            M[i,:] = [3,w]
        elif is_false_negative(a,b,G,Ref):
            M[i,:] = [4,w]
            
    return M 

def count_group_above(data,gid,thresh=-10):
    _data = data[data[:,0] == gid]
    _data = _data[_data[:,1] > thresh]
    return len(_data)

def count_group_below(data,gid,thresh=-10):
    _data = data[data[:,0] == gid]
    _data = _data[_data[:,1] <= thresh]
    return len(_data)


def compute_specificity(m,x,dx=0.1):
    data = np.zeros(len(x))
    for (i,_x) in enumerate(x):
        TN = count_group_below(m,2,thresh=_x)
        FP = count_group_below(m,3,thresh=_x)
        data[i] = 1
        if TN + FP > 0:
            data[i] = float(TN) / (TN + FP)
        #else:
        #    print(x[i])
    return data

def compute_sensitivity(m,x,dx=0.1):
    data = np.zeros(len(x))
    for (i,_x) in enumerate(x):
        TP = count_group_above(m,1,thresh=_x)
        FN = count_group_above(m,4,thresh=_x)
        data[i] = float(TP) / (TP + FN)
    return data

def compute_precision(m,x,dx=0.1):
    data = np.zeros(len(x))
    for (i,_x) in enumerate(x):
        TP = count_group_above(m,1,thresh=_x)
        FP = count_group_above(m,3,thresh=_x)
        if TP + FP > 0:
            data[i] = float(TP) / (TP + FP)
        else: 
            data[i] = 1
    return data


def count_conserved_degrees(deltas,max_deg=4):
    count = np.zeros(max_deg)
    n = float(np.where(deltas>0)[0].size)
    for i in range(1,max_deg+1):
        idx = i - 1
        count[idx] = np.where(deltas == i)[0].size
    return count 


