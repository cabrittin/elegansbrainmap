"""
functions for formatinc igraph data structures

"""

from igraph import VertexClustering
import louvain

def community_vertex_order(community):
    # communinity: Vertex cluster 
    vertex_order = []
    for v in community: vertex_order += v
    return vertex_order


def enforce_community_class(vc,nclass,nid):
    graph = vc.graph
    for (key,vals) in nclass.items():
        comm = set([vc.membership[nid[v]] for v in vals if v in nid])
        if len(comm) < 2: continue
        mod = 0
        for c in comm:
            _membership = vc.membership[:]
            for v in vals: _membership[nid[v]] = c
            _vc = VertexClustering(graph,membership=_membership)
            _mod = _vc.modularity
            if _mod > mod:
                vc = _vc
                mod = vc.modularity
    vc.recalculate_modularity()
                
def get_communities(G,mode=1):
    if mode == 2:
        print('Infomap')
        vc = G.community_infomap(edge_weights='weight')
    elif mode == 3:
        print('Louvain Surprise')
        vc = louvain.find_partition(G,louvain.ModularityVertexPartition)
    elif mode == 4:
        print('Multilevel')
        vc = G.community_multilevel(weights='weight')
    else:
        print('Newman leading eigenvector')
        vc = G.community_leading_eigenvector(weights='weight')

    return vc

def get_network_communities(C,nclass,mode=1):
    nid = dict([(v['name'],v.index) for v in C.A.vs])

    vc_c = get_communities(C.C,mode=mode)
    enforce_community_class(vc_c,nclass,nid)

    vc_e = get_communities(C.E,mode=mode)
    enforce_community_class(vc_e,nclass,nid)

    vc_a = get_communities(C.A,mode=mode)
    enforce_community_class(vc_a,nclass,nid)

    vc_d = get_communities(C.D,mode=mode)
    enforce_community_class(vc_d,nclass,nid)

    return {'vc_c':vc_c,'vc_e':vc_e,'vc_a':vc_a,'vc_d':vc_d}
