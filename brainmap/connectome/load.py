"""
conectome.load.py

Module for routine loading of connectome data

Functions:
---------
from_db(db,chemical=False,electrical=False,adjacency=False,add_poly=False,
            touch_density=False,remove=None,group=None,dataType='igraph'):
   Loads connectome data from database. Returns Connectome object.

scrube_neurons(_neurons)
   Removes neurons in the SCREEN (global variable) list.

created: Christopher Brittin
date: 01 November 2018

"""

import db as db 
from connectome.connectome_igraph import Connectome
from connectome.connectome_networkx import Connectome as nxConnectome
from connectome.connectome_display import Connectome as dConnectome
import networkx as nx

SCREEN = ['old','duplicate','Frag','error','unk']

def from_db(_db,chemical=False,electrical=False,adjacency=False,add_poly=False,
            touch_density=False,remove=None,group=None,dataType='igraph',spatial_domain=0):
    """
    Loads connectome data from database. Returns Connectome object.
    
    Parameters:
    -----------
    _db : str
      database name
    chemical : bool (default False)
      If true, load chemical graph 
    electrical : bool (default False)
      If true, load gap junction graph.
    adjacency : bool (default False)
      If true, load adjacency graph.
    add_poly : bool (default False)
      If true, record number of polyad and monad syanpses.
    touch_density : bool (default (False)
      If true, compute the touch density between adjacent cell in the 
      adjacency graph
    remove : list (default None)
      Remove vertices/cells from the connectome object in the remove list
    group : dict (default None)
      Dictionary (key,val) = (cellName,groupName) to groups cells
    dataType : str (default 'igraph')
      Either use the 'igraph', 'networkx' or 'display' graph data structures. 
    spatial_domain : int
      Restrict adjacency to specified spatial domain: see below

    Returns:
    --------
    Connectome object
    
    """
    
    if _db == 'N2U':
        end = 325
    else:
        end = 500

    con = db.connect.default(_db)
    cur = con.cursor()
    if adjacency:
        neurons = sorted(db.mine.get_adjacency_cells(cur))
        adjacency = db.mine.get_adjacency_data(cur)
        if spatial_domain == 1:
            if _db == 'JSH':
                img = [i.replace('JSHJSH','JSH') for i in db.mine.get_img_number(cur,start=2,end=282)]
            else:
                img = [i.replace('VC','_VC_') for i in db.mine.get_img_number(cur,start=2,end=221)]
                img = [i.replace('NR','_') for i in img]
            adjacency = [a for a in adjacency if a[3] in img]
        elif spatial_domain == 2:
            if _db == 'JSH':
                img = [i.replace('JSHJSH','JSH') for i in db.mine.get_img_number(cur,start=282,end=end)]
            else:
                img = [i.replace('VC','_VC_') for i in db.mine.get_img_number(cur,start=221,end=end)]
                img = [i.replace('NR','_') for i in img]
            adjacency = [a for a in adjacency if a[3] in img]
        elif spatial_domain == 3:
            if _db == 'JSH':
                img = [i.replace('JSHJSH','JSH') for i in db.mine.get_img_number(cur,start=2,end=200)]
            else:
                img = [i.replace('VC','_VC_') for i in db.mine.get_img_number(cur,start=2,end=181)]
                img = [i.replace('NR','_') for i in img]
            adjacency = [a for a in adjacency if a[3] in img]
            
        if dataType == 'igraph':
            C = Connectome(_db,neurons)
        elif dataType == 'networkx':
            C = nxConnectome(_db,neurons)
        elif dataType == 'display':
            C = dConnectome(_db,neurons)
        C.load_adjacency(adjacency,directed=False)
    elif touch_density:
        neurons = sorted(_db.mine.get_adjacency_cells(cur))
        adjacency = db.mine.get_touch_density(cur)
        if dataType == 'igraph':
            C = Connectome(_db,neurons)
        C.load_adjacency(adjacency,directed=True)
    else:
        neurons = sorted(scrub_neurons(db.mine.get_neurons(cur)))
        if dataType == 'igraph':
            C = Connectome(_db,neurons)
        elif dataType == 'networkx':
            C = nxConnectome(_db,neurons)
        elif dataType == 'display':
            C = dConnectome(_db,neurons)
        #C = Connectome(db,neurons)

    if chemical:
        synapses = db.mine.get_synapse_data(cur,'chemical',end=end)
        C.load_chemical(synapses,add_poly=add_poly)

    if electrical:
        synapses = db.mine.get_synapse_data(cur,'electrical',end=end)
        C.load_electrical(synapses)

    if remove: C.remove_cells(remove)
    if group: C.group_cells(group['map'],key=group['key'])
        
    return C


def load_lite(_db,chemical=False,electrical=False,add_poly=False,
              remove=None,group=None,dataType='igraph'):

    """
    Loads connectome data from database. Used to load non-volumetric data.
    Returns Connectome object.
    
    Parameters:
    -----------
    _db : str
      database name
    chemical : bool (default False)
      If true, load chemical graph 
    electrical : bool (default False)
      If true, load gap junction graph.
    add_poly : bool (default False)
      If true, record number of polyad and monad syanpses.
    remove : list (default None)
      Remove vertices/cells from the connectome object in the remove list
    group : dict (default None)
      Dictionary (key,val) = (cellName,groupName) to groups cells
    dataType : str (default 'igraph')
      Either use the 'igraph' or 'networkx' graph data structures. 


    Returns:
    --------
    Connectome object
    
    """    
    con = db.connect.default(_db)
    cur = con.cursor()
    neurons = sorted(scrub_neurons(db.mine.get_neurons(cur)))
    if dataType == 'igraph':
        C = Connectome(_db,neurons)
    elif dataType == 'networkx':
        C = nxConnectome(_db,neurons)
    if chemical:
        synapses = db.mine.get_synapse_data(cur,'chemical')
        C.load_chemical(synapses,add_poly=add_poly)

    if electrical:
        synapses = db.mine.get_synapse_data(cur,'electrical')
        C.load_electrical(synapses)

    if remove: C.remove_cells(remove)
    if group: C.group_cells(group['map'],key=group['key'])
        
    return C    

def filter_synapses(_db,sfilter,args=None,remove=None,group=None,add_poly=False):
    """
    Loads connectome data from database. Returns Connectome object.
    
    Parameters:
    -----------
    _db : str
      database name
    remove : list (default None)
      Remove vertices/cells from the connectome object in the remove list
    group : dict (default None)
      Dictionary (key,val) = (cellName,groupName) to groups cells
    add_poly : bool (default False)
      If true, record number of polyad and monad syanpses.
 
    Returns:
    --------
    Connectome object
    
    """
    
    if _db == 'N2U':
        end = 325
    else:
        end = 500

    con = db.connect.default(_db)
    cur = con.cursor()
    neurons = sorted(scrub_neurons(db.mine.get_neurons(cur)))
    C = nxConnectome(_db,neurons)

    synapses = db.mine.get_synapse_data(cur,'chemical',end=end)
    synapses = sfilter(synapses,args=args)
    #print(synapses)
    C.load_chemical(synapses,add_poly=add_poly)

    if remove: C.remove_cells(remove)
    if group: C.group_cells(group['map'],key=group['key'])
        
    return C

   
def scrub_neurons(_neurons):
    """
    Removes neurons in the SCREEN (global variable) list.

    Parameters:
    ----------
    _neurons : list
     List of cell names to screen

    Returns:
    --------
    List of screened cells
    
    """
    neurons = []
    for n in _neurons:
        remove = False 
        for s in SCREEN:
            if s in n:
                remove = True
                break
        if not remove: neurons.append(n)
    return neurons

def reference_graphs(cfg,min_deg=1,max_deg=4,cfg_section='refgraphs',adj='adj',chem='chem',gap='gap'):
    """
    Loads the reference graphs

    Parameters:
    -----------
    cfg : Configparser, holds the config info
    min_deg: int, default = 1, minimum degree of reproducibility
    max_deg: int, default = 4, maximum degree of reproducibility
    cfg_section: str, default = 'refgraphs', section in config file
    adj: str, default = 'adj', property for adjacency graph
    chem: str, default = 'chem', property for chemical graph
    gap: str, default = 'gap', property for gap junction graph

    Returns:
    ---------
    A : dictionary of reference adjacency graphs by degree (networkx)
    C : dictionary of reference chemical synapses graphs by degree (networkx)
    E : dictionary of reference gap juction graphs by degree (networkx)
    """

    A,C,E = {},{},{}
    for i in range(1,5):
        A[i] = nx.read_graphml(cfg[cfg_section][adj]%i)
        C[i] = nx.read_graphml(cfg[cfg_section][chem]%i)
        E[i] = nx.read_graphml(cfg[cfg_section][gap]%i)

    return A,C,E


def from_graphml(fptr,chemical=False,electrical=False,adjacency=False,remove=None,group=None,**kwargs):
    """
    Loads connectome data from graphml files. Returns Connectome object.
    
    Parameters:
    -----------
    fptr : dict
      Dictionary with file names, format: {'membrane':'path','chem':'path','gap':path}:
    chemical : bool (default False)
      If true, load chemical graph 
    electrical : bool (default False)
      If true, load gap junction graph.
    adjacency : bool (default False)
      If true, load adjacency graph.
    remove : list (default None)
      Remove vertices/cells from the connectome object in the remove list
    group : dict (default None)
      Dictionary (key,val) = (cellName,groupName) to groups cells

    Returns:
    --------
    Connectome object
    
    """
    
    neurons = sorted(scrub_neurons(db.mine.get_neurons(cur)))
    A = nx.read_graphml(fptr['membrane'])
    neurons = A.nodes()

    C = Connectome(cfg['membrane'],neurons)

    if adjacency: C.A = A 
    if chemical: C.C = nx.read_graphml(fptr['chemical'])
    if electrical: C.E = nx.read_graphml(fptr['gap'])

    if remove: C.remove_cells(remove)
    if group: C.group_cells(group['map'],key=group['key'])
        
    return C

