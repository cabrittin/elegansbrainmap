"""
synspecificity.py

Submodule for analyzing synaptic specificity.

Author: Christopher Brittin
Created: 07 February 2018

"""

from scipy.misc import comb
import numpy as np
from scipy.stats import pearsonr
from tabulate import tabulate

import ioaux
from connectome.subcellular import Subcell

def get_bilateral_specificity(C,lrd,left):
    """
    Returns the bilateral (contralateral) specificity probabilities
    
    Parameters:
    -----------
    C : Connecomte object
    lrd : dict
     Left/right dictionary
    left : list
      List of the left nodes

    Returns:
    spec : dict
     Keys are cell names and values are lists with format 
     [gap specificity probability, presynaptic specificity probability,
      postsynaptic specificity probability]
    """
    _pre = bilateral_specificity(C,left,lrd)
    _post = bilateral_specificity(C,left,lrd,mode='post')
    _gap = bilateral_specificity(C,left,lrd,mode='gap') 

    spec = {}
    for nl in left:
        spec[nl] = [_gap[nl].p,_pre[nl].p,_post[nl].p]
    return spec

def get_developmental_specificity(C1,C2,both_nodes):
    """
    Returns the developmental (homologous) specificity probabilities
    
    Parameters:
    -----------
    C1 : Connectome object
    C2 : Connectome
    both_nodes : list
     List of cell names

    Returns:
    spec : dict
     Keys are cell names and values are lists with format 
     [gap specificity probability, presynaptic specificity probability,
      postsynaptic specificity probability]    
    """
    _pre = developmental_specificity(C1,C2,both_nodes)
    _post = developmental_specificity(C1,C2,both_nodes,mode='post')
    _gap = developmental_specificity(C1,C2,both_nodes,mode='gap')

    spec = {}
    for nl in both_nodes:
        spec[nl] = [_gap[nl].p,_pre[nl].p,_post[nl].p]
    return spec

def bilateral_specificity(C,left,lrd,mode='pre'):
    """
    Returns the bilateral (contralateral) specificity probabilities
    for the specified connection type
    
    Parameters:
    -----------
    C : Connecomte object
    lrd : dict
     Left/right dictionary
    left : list
      List of the left nodes
    mode : str
     Connection type ('pre','post','gap'). Default is 'pre'.
    Returns:
    prob : dict
     Keys are cell names and values are Probability objects
    """
    prob = {}
    for nl in left:
        nr = lrd[nl]
        AL = set(C.A.neighbors(nl))
        AR = set([lrd[n] for n in C.A.neighbors(nr)])
        cA = AL & AR
        if mode == 'pre':
            CL = set([n for n in C.C.neighbors(nl) if n in cA])
            CR = set([lrd[n] for n in C.C.neighbors(nr) if lrd[n] in cA])
        elif mode == 'post':
            CL = set([n for n in C.C.predecessors(nl) if n in cA])
            CR = set([lrd[n] for n in C.C.predecessors(nr) if lrd[n] in cA])
        elif mode == 'gap':
            CL,CR = set([]),set([])
            if C.E.has_node(nl):
                CL = set([n for n in C.E.neighbors(nl) if n in cA])
            if C.E.has_node(nr):
                CR = set([lrd[n] for n in C.E.neighbors(nr) if lrd[n] in cA])    
        p = compute_specificity(cA,CL,CR)
        prob[nl] = p
    return prob

def developmental_specificity(C1,C2,nodes,mode='pre'):
    """
    Returns the developmenta (homologous) specificity probabilities
    for the specified connection type
    
    Parameters:
    -----------
    C1 : Connecomte object
    C2 : Connectome object
    nodes : list
     List of cell names
    mode : str
     Connection type ('pre','post','gap'). Default is 'pre'.
    Returns:
    prob : dict
     Keys are cell names and values are Probability objects
    """
    prob = {}
    for _n in sorted(nodes):
        if not(C1.C.has_node(_n) or C2.C.has_node(_n)): continue
        A1 = set(C1.A.neighbors(_n))
        A2 = set(C2.A.neighbors(_n))
        cA = A1 & A2
        S1,S2 = set([]),set([])
        if mode == 'pre':
            S1 = set([n for n in C1.C.neighbors(_n) if n in cA])
            S2 = set([n for n in C2.C.neighbors(_n) if n in cA])
        elif mode == 'post':
            S1 = set([n for n in C1.C.predecessors(_n) if n in cA])
            S2 = set([n for n in C2.C.predecessors(_n) if n in cA])
        elif mode == 'gap':
            if C1.E.has_node(_n):
                S1 = set([n for n in C1.E.neighbors(_n) if n in cA])
            if C2.E.has_node(_n):
                S2 = set([n for n in C2.E.neighbors(_n) if n in cA])
        p = compute_specificity(cA,S1,S2)
        prob[_n] = p
    return prob
        
def compute_specificity(A,C1,C2):
    """
    Computes the specificty probability

    Parameters:
    -----------
    A : list
     List of common (physical) neighbors
    C1 : list
     List of synaptic partners in A
    C2 : list
     List of synaptic partners in A

    Returns:
     Probability object
    """
    M = len(A)
    cC = C1 & C2
    k = len(cC)
    K = min(len(C1),len(C2))
    tmp1 = len(C1)
    tmp2 = len(C2)
    c1 = max(tmp1,tmp2)
    c2 = min(tmp1,tmp2)
    m = M - k - c1 - c2
    den1 = comb(M,c1)
    den2 = comb(M,c2)
    den = den1*den2
    
    num = 0

    for _k in range(k,K+1):
        _c1 = c1 - _k
        _c2 = c2 - _k
        _m = M - _k - _c1 - _c2
        
        num1 = comb(M,_k)
        num2 = comb(M-_k,_c1)
        num3 = comb(M-_k-_c1,_c2)
        #num4 = comb(M-_k-_c1-_c2,_m)

        #print(M,_k+_c1+_c2+_m,_k,_c1,_c2,_m,num1,num2,num3,num4)
        num += num1*num2*num3#*num4
    return Prob(M,k,c1,c2,m,num / den)


def synapse_positions(cur,neuron):
    """
    Returns synapse positions

    Parameters:
    -----------
    cur : MySQLdb cursor
    neuron : str
      cell name

    Returns:
    S : dict
     Three entries 'pre', 'post' and 'gap', which have the synapse positions
     for each synapse type. 
    """
    S = {'pre':{},'post':{},'gap':{}}

    get_synapses_pos(cur,neuron,S['pre'],'pre','post','chemical')
    get_synapses_pos(cur,neuron,S['post'],'post','pre','chemical')
    get_synapses_pos(cur,neuron,S['gap'],'pre','post','electrical')
    get_synapses_pos(cur,neuron,S['gap'],'post','pre','electrical')

    for k in S:compute_average_pos(S[k])

    return S

def format_left_right_subcell(S1,S2,lrd):
    """
    Formats the data structures for bilateral (contralateral)
    comaprisons
    
    Parameters:
    -----------
    S1 : dict
     Dictionary of synapses posisions
    S2 : dict
     Dictionary of synapse positions
    lrd : dict
     Left/right dictionary
    """
    S1['pre'],S2['pre'] = format_left_right_dicts(S1['pre'],S2['pre'],lrd)
    S1['post'],S2['post'] = format_left_right_dicts(S1['post'],S2['post'],lrd)
    S1['gap'],S2['gap'] = format_left_right_dicts(S1['gap'],S2['gap'],lrd)

    return S1,S2

def format_left_right_dicts(l1,l2,lrd):
    _l2new = dict([(lrd[l],l2[l]) for l in l2 if l in lrd])

    common = set(l1.keys()) & set(_l2new.keys())
    l1new,l2new = {},{}
    for c in common:
        l1new[c] = l1[c]
        l2new[c] = _l2new[c]
    if 'mean' in l1: l1new['mean'] = l1['mean']
    if 'mean' in l2: l2new['mean'] = l2['mean']
    if 'std' in l1: l1new['std'] = l1['std']
    if 'std' in l2: l2new['std'] = l2['std']
    if 'size' in l1: l1new['size'] = l1['size']
    if 'size' in l2: l2new['size'] = l2['size']
    return l1new,l2new

def format_developmental_subcell(S1,S2):
    S1['pre'],S2['pre'] = format_developmental_dicts(S1['pre'],S2['pre'])
    S1['post'],S2['post'] = format_developmental_dicts(S1['post'],S2['post'])
    S1['gap'],S2['gap'] = format_developmental_dicts(S1['gap'],S2['gap'])    

    return S1,S2
    
def format_developmental_dicts(d1,d2):
    common = set(d1.keys()) & set(d2.keys())
    d1new,d2new = {},{}
    for c in common:
        d1new[c] = d1[c]
        d2new[c] = d2[c]
    if 'mean' in d1: d1new['mean'] = d1['mean']
    if 'mean' in d2: d2new['mean'] = d2['mean']
    return d1new,d2new

def make_subcell_table(S1,S2):
    """
    Will display the synapse positions in table format
    """
    def add_data(col,S1,S2,data):
        common = set(S1[col]) | set(S2[col])
        for c in common:
            a,b = 'NA','NA'
            if c in S1[col]: a = S1[col][c]
            if c in S2[col]: b = S2[col][c]
            data.append([col,c,a,b])

    
    data = []
    add_data('pre',S1,S2,data)
    add_data('post',S1,S2,data)
    add_data('gap',S1,S2,data)
    print(tabulate(data,
                   headers=['Mode','Left pos','Right pos'],
                   tablefmt='orgtbl'))

def get_synapses_pos(cur,neuron,pos,fromCol,toCol,stype):
    """
    Updates synapses positions for neurons
    Decided not to place in the db module because
    the function returns a more sophiticated data type.

    Parameters:
    -----------
    cur : MySQLdb cursor
    neuron : str
     cell name
    pos : dict
     Synapses positons are written here. pos will have the 
     forma pos[syn_parnter] = {'loc':[], 'weight':[]}, where
     syn_partner is the name of the synaptic partner, loc is the 
     location of all synapses, weight is the lenght of the synapse 
     in EM sections. 
    fromCol : str
     Specifies the presynaptic column. Used to accomodate gap junctions
     which are undirected and do not have pre/post partners
    toCol : str
      Specified the postsynaptic column. Used to accomodate gap junctions
     which are undirected and do not have pre/post partners
    stype : str
      Synapse type ('chemical', 'electrical')
    
    Note that 'pos' variable is update, but no value is returned. 

    """

    sql = ("select %s,sections,position,continNum "
           "from synapsecombined "
           "join synapseskeleton "
           "on synapsecombined.mid = synapseskeleton.mid "
           "where synapseskeleton.neuron = '%s' "
           "and synapsecombined.type = '%s' "
           "and synapsecombined.%s like '%%%s%%'"
           %(toCol,neuron,stype,fromCol,neuron))
    cur.execute(sql)
    for a in cur.fetchall():
        nodes = a[0].split(',')
        for n in nodes:
            if n not in pos: pos[n] = {'loc':[],'weight':[]}
            pos[n]['loc'].append(a[2])
            pos[n]['weight'].append(a[1])

    
def compute_average_pos(syn):
    """
    Computes the average position of the synapses
    
    Parameters:
    -----------
    syn : dict
     In the formatted data structure. Keys are cell names. 
     Values are dicts(['loc':[...],'weights':[...]])
    
    Note that syn is updated but not returned. 
    Will add the (key, values):
     'mean' : weighted mean syanpse position
     'std'  : weighted standard deviation
     'size' : number of synapses
    
    """
    all_syn_loc = []
    all_syn_weight = []
    for s in syn:
        all_syn_loc += syn[s]['loc']
        all_syn_weight += syn[s]['weight']
        ###Updates list of synapses with the weighted average
        syn[s] = np.average(syn[s]['loc'],weights=syn[s]['weight'])

    if all_syn_weight:
        mu = np.average(all_syn_loc,weights=all_syn_weight)
        syn['mean'] = mu
        loc = np.array(all_syn_loc)
        variance = np.average((loc-mu)**2,weights=all_syn_weight)
        std = np.sqrt(variance)
        #w = np.array(all_syn_weight)
        #num = np.sum(w*(loc-syn['mean'])**2)
        N = float(len(all_syn_weight))
        #if N > 1:
        #    den = (N-1)/N*(np.sum(w))
        #else:
        #    den = np.sum(w)
        syn['std'] = std
        syn['size'] = N
    else:
        syn['mean'] =-999
        syn['std'] = -999
        syn['size'] = -999
        
def get_outliers(pre,post,ndict):
    prob_both = []
    outliers = []
    count1,count2,count3,count4 = 0,0,0,0
    for n in sorted(pre):
        prob_both.append((pre[n].p,post[n].p))
        if (pre[n].p > 0.05 and post[n].p > 0.05):
            count4 += 1
            print(n,
                  pre[n].p,post[n].p,
                  pre[n].c1,post[n].c1,
                  pre[n].k,post[n].k,
                  pre[n].M,post[n].M)
            outliers.append([ndict[n],pre[n].p,post[n].p])
        elif pre[n].p <= 0.05 and post[n].p <= 0.05:
            count1 += 1
        elif pre[n].p <= 0.05 and post[n].p > 0.05:
            count2 += 1
        elif pre[n].p > 0.05 and post[n].p <= 0.05:
            count3 += 1
    N = float(len(pre))
    print('Type1: %1.4f, Type2: %1.4f, Type3: %1.4f, Type4: %1.4f' %(count1/N,count2/N,count3/N,count4/N))
    print('Type1: %d, Type2: %d, Type3: %d, Type4: %d, Total: %d' %(count1,count2,count3,count4,N))
    return prob_both,outliers,[count1/N,count2/N,count3/N,count4/N]

def plot_hist(ax,data,bins,col,label):
    ax.hist(data,bins=bins,range=(0,1),
            histtype='step',linewidth=4,
            cumulative=1,
            normed=1,
            color=col,
            label = label)    

def write_specificity(fout,pre,post,ndict):
    d = {}
    for n in pre:
        d[ndict[n]] = [pre[n].p,post[n].p]

    ioaux.write.from_dict(fout,d)
    
    
class Prob:
    def __init__(self,M,k,c1,c2,m,p):
        self.M = M
        self.k = k
        self.c1 = c1
        self.c2 = c2
        self.m = m
        self.p = p

    def get_attrib(self):
        return [self.M,self.k,self.c1,self.c2,self.p]

def get_source_data(cur,cells):
    data = {}
    for n in cells:
        data[n] = synapse_positions(cur,n)
    return data
            
def get_bilateral_subcell_specificity(cur,neurons,lrd):
    """
    Returns the bilateral subcellular specificity values

    Parameters:
    cur : MySQLdb cursor
    neurons : str
     list of cell names
    lrd : dict
     Left/right dictionary

    Returns:
    --------
    delta : dict
     (key,value) = (cell, [gap_delta,pre_delta,post_delta])
    """
    delta = {}
    for [nl,nr] in neurons:
        #print(nl)
        Sl = synapse_positions(cur,nl)
        Sr = synapse_positions(cur,nr)
        Sl,Sr = format_left_right_subcell(Sl,Sr,lrd)
        delta[nl] = [[],[],[]]
        #print('\tpre')
        delta[nl][1] = compute_delta(Sl['pre'],Sr['pre'])
        #print('\tpost')
        delta[nl][2] = compute_delta(Sl['post'],Sr['post'])
        #print('\tgap')
        delta[nl][0] = compute_delta(Sl['gap'],Sr['gap'])

    return delta

def get_developmental_subcell_specificity(cur1,cur2,both_nodes=None):
    """
    Returns the developmental subcellular specificity values

    Parameters:
    cur1 : MySQLdb cursor
    cur2 : MySQLdb cursor
    both_nodes : str
     list of cell names

    Returns:
    --------
    delta : dict
     (key,value) = (cell, [gap_delta,pre_delta,post_delta])
    """    

    delta = {}
    for n in sorted(both_nodes):
        #print(n)
        S1 = synapse_positions(cur1,n)
        S2 = synapse_positions(cur2,n)
        S1,S2 = format_developmental_subcell(S1,S2)
        delta[n] = [[],[],[]]
        #print('\tpre dev')
        delta[n][1] = compute_delta(S1['pre'],S2['pre'])
        #print('\tpost dev')
        delta[n][2] = compute_delta(S1['post'],S2['post'])
        #print('\tgap dev')
        delta[n][0] = compute_delta(S1['gap'],S2['gap'])

    return delta    
    
    
def batch_display(corr):
    data = []
    for n in sorted(corr):
        data.append( [n] + corr[n])
    print(tabulate(data,
                   headers=['Neuron','Presyn','Postsyn','Gap'],
                   tablefmt='orgtbl'))    


def compute_delta(l1,l2):
    """
    Returns list of differences between the mean syanpse posistions
    of synaptic partners
    
    Parameters:
    l1 : dict
     Dictionary with (key,value) = (synapse partner, average position)
    l2 : dict 
     Dictionary with (key,value) = (synapse partner, average position)
    
    Returns:
    delta : list
     list of difference in synapses positons

    Note that keys: 'mean', 'std', 'size' are skipped
    """
    delta = []
    mean_delta = []
    N1,N2 = l1['size'],l2['size']
    sig1,sig2 = l1['std'],l2['std']
    #compute average standard deviation
    if (N1 > 1 and N2 > 1):
        num = (N1-1)*sig1**2 + (N2-1)*sig2**2
        den = (N1 + N2 -2)
    else:
        num = sig1**2 + sig2**2
        den = 2.
    
    sp = np.sqrt(float(num)/den)
    for n in l1:
        if n == 'mean': continue
        if n == 'std': continue
        if n == 'size': continue
        _delta = abs(l1[n] - l2[n])
        #print('\tstd: %s (%1.2f,%1.2f), %1.2f'%(n,l1[n],l2[n],sp))
        #if _delta > 0.5:
        #    print('\tcheck: %s (%1.2f,%1.2f)'%(n,l1[n],l2[n]))
        delta.append((l1[n] - l2[n]))
        #mean_delta.append(abs(l1[n] - l2['mean']))
    return delta

def compute_diff(l1,l2):
    return np.mean([l1[n]-l2[n] for n in l1])



    
