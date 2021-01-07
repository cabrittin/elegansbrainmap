""""
random.py

Modules for generating randomized networks.

Functions
---------
edge_switch(G,iters=1000)
  Randomly switch edges in graph G (igraph). Do iters random switches.
  Preserves in degree and out degree.

mp_edge_switch(G,proFunc,njobs=2,sub_iters=500,es_iters=1000):
  Do random edge switches in parallel

_mp_edge_switch(G,procFunc,output,sub_iters,es_iters):
  Target function for mp_edge_switch.


"""
from random import randint
from random import sample

import multiprocessing as mp
import progressbar


def edge_switch(G,iters=1000):
    """
    Randomly switch edges in graph G (igraph). Do iters random switches.
    Preserves in degree and out degree.

    Parameters
    ----------
    G : iGraph
     Graph in which edges will be randomly switched
    iters : int (default 1000)
     Number of edge switches to perform 

    """
    N = G.ecount() - 2
    idx = 0
    while idx < iters:
        i = randint(0,N)
        j = i
        while j == i: j = randint(0,N)
        ei = G.es[i]
        ej = G.es[j]
        u1 = ei.source
        v1 = ei.target
        u2 = ej.source
        v2 = ej.target
        iattr = ei.attributes()
        jattr = ej.attributes()
        cond1 = len(set([u1,v1,u2,v2])) == 4
        cond2 = not G.are_connected(u1,v2) 
        cond3 = not G.are_connected(u2,v1) 
        if cond1 and cond2 and cond3: 
            G.delete_edges([i,j])
            G.add_edge(u1,v2,**iattr)
            G.add_edge(u2,v1,**jattr)
            idx += 1
        
    
def mp_edge_switch(G,procFunc,njobs=2,sub_iters=500,es_iters=1000):
    """
    Do random edge switches in parallel
    
    Parameters
    ----------
    G : iGraph
     Graph object
    procFunc : function
     User supplied function for how to handle the results
    njobs : int
     Number cpus to use
    sub_iters : int
     number of random graph to generate per job
    es_iters : int
     Number of edge switches in a single graph
    """
    output = mp.Queue()
    processes = [mp.Process(target=_mp_edge_switch,
                            args=(G,procFunc,output,sub_iters,es_iters))
                 for i in range(njobs)]
    
    for p in processes:
        p.start()

    for p in processes:
        p.join()

    return [output.get() for p in processes]


def _mp_edge_switch(G,procFunc,output,sub_iters,es_iters):
    """
    Target function for mp_edge_switch
    
    Parameters
    ----------
    G : iGraph
     Graph object
    procFunc : function
     User supplied function for how to handle the results
    output : multprocessor.Queue
     Queue to output results
    sub_iters : int
     number of random graph to generate per job
    es_iters : int
     Number of edge switches in a single graph
    """   

    bar = progressbar.ProgressBar(maxval=sub_iters,
                                  widgets=[progressbar.Bar('.','[',']'),
                                           'Edge switch randomization',
                                           progressbar.Percentage()])
    bar.start()
    R = G.copy()
    results = []
    for i in range(sub_iters):
        edge_switch(R,iters=es_iters)
        results.append(procFunc(R))
        bar.update(i+1)
    output.put(results)
