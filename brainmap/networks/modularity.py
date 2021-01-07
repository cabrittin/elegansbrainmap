from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np

def modularity_info_metric(vc,metric):
    ca = metric(vc['vc_c'],vc['vc_a'])
    ea = metric(vc['vc_e'],vc['vc_a'])
    da = metric(vc['vc_d'],vc['vc_a'])

    return [ca,ea,da]

def modularity_quality_metric(vc,metric,
                              key_order=['vc_c','vc_e','vc_d','vc_a']):
    return [metric(vc[k]) for k in key_order]

########################
#  Information scores  #
########################

def normalized_mutual_info(comm1,comm2):
    """
    comm1: community 1
    comm2: community 2
    
    list of community labels, index i is integer community label of node i 
    """

    return normalized_mutual_info_score(comm1.membership,comm2.membership)

def adjusted_mutual_info(comm1,comm2):
    """
    comm1: community 1
    comm2: community 2
    
    list of community labels, index i is integer community label of node i 
    """

    return adjusted_mutual_info_score(comm1.membership,comm2.membership)


def adjusted_rand(comm1,comm2):
    """
    comm1: community 1
    comm2: community 2
    
    list of community labels, index i is integer community label of node i 
    """

    return adjusted_rand_score(comm1.membership,comm2.membership)


########################
#  Quality scores      #
########################

def modularity(comm):
    """
    comm: community 
    """
    return comm.modularity
    
    
def coverage(comm,edge_attr=None):
    """
    comm: community
    """
    num,den = 0.,0.
    for e in comm.graph.es:
        w = 1
        if edge_attr: w = e[edge_attr]
        den += w
        i = e.source
        j = e.target
        if comm.membership[i] == comm.membership[j]:
            num += w
    return num/den            

def conductance(comm):
    """
    comm: community
    """
    A = comm.graph.get_numpy_array(edge_attr='weight')
    k = list(set(comm.membership))
    membership = np.array(comm.membership)
    phi = np.zeros(len(k))
    for _k in k:
        sk = np.where(membership ==_k)[0]
        nsk = np.where(membership !=_k)[0]
        num = np.sum(A[sk][:,nsk])
        denA = np.sum(A[sk,:]) - np.sum(A[sk,sk])
        denB = np.sum(A[nsk,:]) - np.sum(A[nsk,nsk])
        den = min(denA,denB)
        if den == 0:
            phi[_k] = 0
        else:
            phi[_k] = num/float(den)
    return 1 - np.mean(phi)
            
