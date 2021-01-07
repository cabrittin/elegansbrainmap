import numpy as np

def arcsine(data):
    return 2*np.arcsin(np.sqrt(data))

def get_cf(C,vertices=None,_arcsine=False):
    if _arcsine:
        data =  {'pre':arcsine(get_cpre(C,vertices=vertices)),
                 'post':arcsine(get_cpost(C,vertices=vertices)),
                 'gap':arcsine(get_cgap(C,vertices=vertices))}
        if C.D:
            data['all'] = arcsine(get_call(C,vertices=vertices))
        
    else:
        data =  {'pre':get_cpre(C,vertices=vertices),
                 'post':get_cpost(C,vertices=vertices),
                 'gap':get_cgap(C,vertices=vertices)}
        if C.D:
           data['all'] = get_call(C,vertices=vertices) 
    return data

def get_call(C,vertices=None):
    if not C.D: return None
    if vertices:
        v = vertices
    else:
        v = [v['name'] for v in C.A.vs]
    con = np.array(C.D.degree(v,mode='ALL',loops=False))
    adj = np.array(C.A.degree(v,mode='ALL',loops=False))
    return con/adj

def get_cpre(C,vertices=None):
    if vertices:
        v = vertices
    else:
        v = [v['name'] for v in C.A.vs]
    con = np.array(C.C.degree(v,mode='OUT',loops=False))
    adj = np.array(C.A.degree(v,mode='ALL',loops=False))
    return con/adj

def get_cpost(C,vertices=None):
    if vertices:
        v = vertices
    else:
        v = [v['name'] for v in C.A.vs]
    con = np.array(C.C.degree(v,mode='IN',loops=False))
    adj = np.array(C.A.degree(v,mode='ALL',loops=False))
    return con/adj

def get_cgap(C,vertices=None):
    if vertices:
        v = vertices
    else:
        v = [v['name'] for v in C.A.vs]
    con = np.array(C.E.degree(v,mode='ALL',loops=False))
    adj = np.array(C.A.degree(v,mode='ALL',loops=False))
    return con/adj

def get_adj_deg(C,vertices=None,norm=False):
    return get_deg(C.A,vertices=vertices,mode='All',norm=norm)

def get_deg(G,vertices=None,mode='All',norm=False):
    if vertices:
        _vertices = [v['name'] for v in G.vs]
        v = [_v for _v in vertices if _v in _vertices]
    else:
        v = [v['name'] for v in G.vs]

    deg = {}
    if norm:
        N = float(G.vcount())
        for _v in v:
            deg[_v] = G.degree(_v,mode=mode,loops=False) / N
    else:
        for _v in v:
            deg[_v] = G.degree(_v,mode=mode,loops=False)
    return deg
    
def get_td_bounds(td):
    td = np.log(td)
    mu = np.mean(td)
    std2 = 2*np.std(td)
    return [mu - std2, mu + std2] 

def get_corresponding_edge_attr(G1,G2,attr1='weight',attr2='weight'):
    _attr1,_attr2 = [],[]
    _edges = []
    for _u in G1.vs:
        u = _u['name']
        neigh1 = G1.vs[G1.neighbors(u,mode='out')]['name']
        neigh2 = G2.vs[G2.neighbors(u,mode='out')]['name']
        neigh = set(neigh1) | set(neigh2)
        u1 = _u.index
        u2 = G2.vs.find(name=u).index
        for v in neigh:
            v1 = G1.vs.find(name=v).index
            v2 = G2.vs.find(name=v).index
            w1,w2 = 0,0
            if G1.are_connected(u1,v1):
                w1 = G1.es[G1.get_eid(u1,v1)][attr1]                
            if G2.are_connected(u2,v2):
                w2 = G2.es[G2.get_eid(u2,v2)][attr2]
            _edges.append(u + '-' + v)
            _attr1.append(w1)
            _attr2.append(w2)
        
    return _edges,_attr1,_attr2

def get_corresponding_out_strength(G1,G2,edge_attr='weight'):
    attr1,attr2 = [],[]
    for _u in G1.vs:
        u = _u['name']
        neigh1 = G1.vs[G1.neighbors(u,mode='out')]['name']
        neigh2 = G2.vs[G2.neighbors(u,mode='out')]['name']
        neigh = set(neigh1) | set(neigh2)
        u1 = _u.index
        u2 = G2.vs.find(name=u).index
        s1 = float(G1.vs[u1]['out_strength'])
        for v in neigh:
            v1 = G1.vs.find(name=v).index
            v2 = G2.vs.find(name=v).index
            w1,w2 = 0,0
            if G1.are_connected(u1,v1):
                w1 = G1.es[G1.get_eid(u1,v1)][attr]                
            if G2.are_connected(u2,v2):
                w2 = G2.es[G2.get_eid(u2,v2)][attr]
            attr1.append(w1)
            attr2.append(w2)
        
    return attr1,attr2     

def filter_corresponding_tds(edges,td1,td2,tdbounds):
    [tmin1,tmax1] = tdbounds[0]
    [tmin2,tmax2] = tdbounds[1]
    N = len(td1)
    _edges,_td1,_td2 =[], [],[]
    
    for i in range(N):
        if td1[i] == 0 or td2[i] == 0: continue
        log1 = np.log(td1[i])
        log2 = np.log(td2[i])
        if ((log1 > tmin1 and log1 < tmax1)
            and (log2 > tmin2 and log2 < tmax2)):
            _edges.append(edges[i])
            _td1.append(log1)
            _td2.append(log2)
    return _edges,_td1,_td2
            
def get_neighborhood_similarity(G1,G2,vertices,mode='out'):
    """
    Gives -1 value for vertices not found in graph
    """
    
    for u in vertices:
        try:
            neigh1 = G1.vs[G1.neighbors(u,mode=mode)]['name']
            neigh2 = G2.vs[G2.neighbors(u,mode=mode)]['name']
            nsim = neighborhood_similarity(set(neigh1),set(neigh2))
        except ValueError:
            nsim = -1
    
        yield (u,nsim)

def get_neighborhood_overlap_similarity(G,vertices,mode='out'):
    """
    Gives -1 value for vertices not found in graph
    """
    for u in vertices:
        try:
            neigh1 = set(G.vs[G.neighbors(u,mode=mode)]['name'])
            sim = -1
            for v in neigh1:
                neigh2 = set(G.vs[G.neighbors(v,mode=mode)]['name'])
                _sim = neighborhood_similarity(neigh1,neigh2)
                sim = max(sim,_sim)
        except ValueError:
            sim = -1

        yield (u,sim)
            
                      
def neighborhood_similarity(neigh1,neigh2):
    #return 1 - jaccard(neigh1,neigh2)
    return jaccard(neigh1,neigh2)
    
def jaccard(A,B):
    """
    Compute Jaccard similarity between sets A and B
    """
    intersect = len(A & B)
    union = len(A | B)
    if union == 0:
        return 0
    else:
        return intersect / float(union)

def compute_similarity_score(score):
    idx = np.where(score[:,0] > 0)[0]
    score = score[idx,:]
    score = 2*score - 1
    score = np.arctan(score)
    std = np.std(score[:,0])
    _score = (score[:,1] - score[:,0])/std
    return _score
    #return np.log(score[:,0] / score[:,1])


def get_adj_poly(C,A):
    """
    Get adjacency weigths labeled as: 
        no synapse => 0 
        polyadic   => 1
        monadic    => 2
    
    Input:
       C: chemical network
       A: adjacency network
    
    Return:
       data: [adj weight, synapse type]

    """
    SCALE = 5*90*(1e-6)
    C.to_undirected(combine_edges=sum)
    N = A.ecount()
    data = np.zeros((N,2))
    for i in range(N):
        e = A.es[i]
        data[i,0] = e['weight']
        if C.are_connected(e.source,e.target):
            sp = C.es[C.get_eid(e.source,e.target)]['Sp']
            if sp > 0:
                data[i,1] = 1
            else:
                data[i,1] = 2
        
    data[:,0] = np.log(data[:,0]*SCALE)
    return data    

def get_adj_poly_data(C):
    data = get_adj_poly(C.C,C.A)
    _mon = np.where(data[:,1] == 2)[0]
    _poly = np.where(data[:,1] == 1)[0]
    _zero = np.where(data[:,1] == 0)[0]
    return [data[_zero,0],data[_poly,0],data[_mon,0]]
    
def eff_diff(data1,data2,pool=True,std_idx = 0):
    if pool:
        data3 = []
        data3.extend(data1)
        data3.extend(data2)
        std = np.std(data3)
    elif std_idx == 0:
        std = np.std(data1)
    else:
        std = np.std(data2)
    diff = np.array(data1) - np.array(data2)
    return diff / std

def get_venn_data(G1,G2):
    g1only = 0
    both = 0
    for (u,v) in G1.edges():
        if G2.has_edge(u,v):
            both += 1
        else:
            g1only += 1
    g2only = G2.number_of_edges() - both
    return [g1only,g2only,both]
