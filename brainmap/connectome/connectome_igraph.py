"""
connectome.py

Connectome data structures. Inherits from iGraph

See http://igraph.org/python/

Required 3rd party packages:
  igraph
  csv
  numpy

Author: Christopher Brittin

"""

import igraph
import csv
import numpy as np

import ioaux

SCREEN = ['old']

class Connectome:
    """
    Class to represent connectome data.
    
    ...

    Attributes
    ----------
    db : str
      database name
    size : int
      number of neurons in graph
    neurons : list
      list of neuron names
    C  : Network
      Chemical connectivity graph
    E  : Network
      Gap junction connectivity graph
    A  : Network
      Adjacency (physical) connectivity graph
    D  : Network
      Combined chemical and gap junction graph

    
    Methods
    --------
    update_cells(_neurons)
      Update the neuron list

    remove_self_loops()
      Remove self loops from graphs C and E

    remove_cells(vertices)
      Remove vertices from graphs C, E and A

    group_cells(groups,key='group')
      Group vertices based on dictionary groups. The grouping identified
      by key (default 'group'). So multiple groups can be assigned to 
      same graphs.
    
    load_chemical(synapses,add_poly=False)
      Create chemical connectivity graph by loading edges from synapses. 
      If add_poly, then the number of polyad and monad synapses is tracked.

    load_electrical(synapses,add_poly=False)
      Create gap junction connectivity graph by loading edges from synapses. 
      If add_poly, then the number of polyad and monad synapses is tracked.

    load_edges(G,vertices,edges,add_poly=False)
      Load edges between vertices into graph G. If add_poly, then the
      number of polyad and monad synapses is tracked.

    load_adjacency(_adjacency,directed=False)
      Load adjacency graph from _adjacency edges. If directed, adjacency graph
      will be made directed. 

    combine_chem_and_elec()
      Combine the chemical and gap junction connectivity graphs

    reduce_to_adjacency()
      Reduce chemical and gap junction connectivity graphs to nodes and
      edges found in the adjacency graph

    

    """
    
    def __init__(self,db,neurons):
        """
        Parameters:
        -----------
        db : str
           database name
        neurons : list
           list of neuron/cell/node names
        """
        self.db = db
        self.size = len(neurons)
        self.neurons = neurons
        self.C = None
        self.E = None
        self.A = None
        self.D = None

    def update_cells(self,_neurons):
        """
        Update the neuron/cell/node list

        Parameters:
        -----------
        _neurons : list 
           list of neurons/cell/node/names

        """
        self.neurons = _neurons
        self.size = len(_neurons)

    def remove_self_loops(self):
        """
        Remove self loops from graph C and E
        """
        if self.C: self.C.simplify(multiple=True,loops=True)
        if self.E: self.E.simplify(multiple=True,loops=True)
        
    def remove_cells(self,vertices):
        """"
        Remove vertices from graphs C, E and A
        
        Parameters
        ----------
        vertices : list
          List of vertex names

        """
        self.neurons = list(set(self.neurons) - set(vertices))
        if self.C:self.C.remove_vertices(vertices)
        if self.E:self.E.remove_vertices(vertices)
        if self.A:self.A.remove_vertices(vertices)

    def group_cells(self,groups,key='group'):
        """
        Group vertices based on list groups. The grouping identified
        by key (default 'group'). So multiple groups can be assigned to 
        same graphs.                
        
        Parameters
        ----------
        groups : list
          List of vertex memberships. Order of list 
          should correspond to index labels of vertices.
        
        key : str
          Group ID. Multiple groups can be assigned and
          retrieved with the ID.

        """
        if self.C:
            self.C.assign_membership(groups,key=key)
            self.C.group_vertices(key)
        if self.E:
            self.E.assign_membership(groups,key=key)
            self.E.group_vertices(key)
        if self.A:
            self.A.assign_membership(groups,key=key)
            self.A.group_vertices(key)
                
    def load_chemical(self,synapses,add_poly=False):
        """
        Create chemical connectivity graph by loading edges from synapses. 
        If add_poly, then the number of polyad and monad synapses is tracked.       
        
        Parameters:
        synapses : list
          list of synapse data with row format
          pre_cell,post_cell(s),synapse_weight,synapse_id,data_series
        add_poly : bool (default False)
          If true, then tabulate polyadic and monadic synapses
        
        """
        self.C = Network(directed=True)
        self.C.add_vertices(self.neurons)
        self.load_edges(self.C,self.neurons,synapses,add_poly=add_poly)

    def load_electrical(self,synapses,add_poly=False):
        """
        Create gap junction connectivity graph by loading edges from synapses. 
        If add_poly, then the number of polyad and monad synapses is tracked.       
        
        Parameters:
        synapses : list
          list of synapse data with row format
          pre_cell,post_cell(s),synapse_weight,synapse_id,data_series
        add_poly : bool (default False)
          If true, then tabulate polyadic and monadic synapses

        """
        self.E = Network()
        self.E.add_vertices(self.neurons)
        self.load_edges(self.E,self.neurons,synapses,add_poly=add_poly)
                
    def load_edges(self,G,vertices,edges,add_poly=False):
        """
        Load edges between vertices into graph G. If add_poly, then the
        number of polyad and monad synapses is tracked. 

        Parameters:
        -----------
        G : Network
          Graph into which edges will be loaded
        vertices : list 
          list of vertex names. At least one vertex in edge
          must be in the list vertex names
        edges : list
          list of edge data with row format
          pre_cell,post_cell(s),synapse_weight,synapse_id,data_series
        add_poly : bool (default False)
          If true, then tabulate polyadic and monadic synapses
        """
        eid = -1
        for e in edges:
            pre = ioaux.format.rm_brack(e[0])
            if pre not in vertices: continue
            #i_pre = self.neurons[pre]
            _post = list(set(map(ioaux.format.rm_brack,e[1].split(','))))
            if add_poly:
                if len(_post) == 1:
                    poly = 'S'
                else:
                    poly = 'Sp'
            if self.db == 'N2U' and e[4] in ['VC','DC']:
                w = int(e[2])#2*int(e[2]) - 1
            else:
                w = int(e[2])
                
            for post in _post:
                if post not in vertices: continue 
                #i_post = self.neurons[post]
                if not G.are_connected(pre,post):
                    eid += 1
                    G.add_edges([(pre,post)])
                    G.es[eid]['weight'] = 0
                    G.es[eid]['count'] = 0
                    if add_poly:
                        G.es[eid]['S'] = 0
                        G.es[eid]['Sp'] = 0
                _eid = G.get_eid(pre,post)
                G.es[_eid]['weight'] += w
                G.es[_eid]['count'] += 1
                if add_poly:
                    G.es[_eid][poly] += 1
        #G.vs.select(_degree=0).delete()

    def load_adjacency(self,adjacency,directed=False):
        """
        Load adjacency graph from _adjacency edges. If directed, adjacency graph
        will be made directed.        

        Parameters:
        -----------
        adjacency : list
          List of adjacency data with row format
          cell1,cell2,amount_of_contact,section_number
        directed: bool
          If true, the adjacency graph will be directed
        
        """
        self.A = Network(directed=directed)
        self.A.add_vertices(self.neurons)
        eid = -1
        for (i,j,weight,imgNum) in adjacency:
            weight = weight
            count = 1
            if self.db == 'N2U' and 'VC' in imgNum:
                weight *=1#2
                count = 1
            if not self.A.are_connected(i,j):
                eid += 1
                self.A.add_edges([(i,j)])
                self.A.es[eid]['weight'] = 0
                self.A.es[eid]['count'] = 0
            _eid = self.A.get_eid(i,j)
            self.A.es[_eid]['weight'] += weight
            self.A.es[_eid]['count'] += count

    def combine_chem_and_elec(self):
        """
        Combine the chemical and gap junction connectivity graphs
        Combined graph stored in attribute self.D

        """
        self.D = self.C.copy()
        for e in self.E.es:
            i = self.E.vs[e.source]['name']
            j = self.E.vs[e.target]['name']
            w = 0.5*e['weight']
            c = 0.5*e['count']
            if self.D.are_connected(i,j):
                eid = self.D.get_eid(i,j)
                self.D.es[eid]['weight'] += w
                self.D.es[eid]['count'] += c
            else:
                self.D.add_edge(i,j,weight=w,count=c)

    def reduce_to_adjacency(self):
        """
        Reduce chemical and gap junction connectivity graphs to nodes and
        edges found in the adjacency graph  
        """
        if self.C: self.C.reduce_to(self.A)
        if self.E: self.E.reduce_to(self.A)
        
class Network(igraph.Graph):
    """
    Class for storing graph data. Inherits class Graph
    from package igraph. See http://igraph.org/python/ 
    for attributes and API.
    
    Methods
    -------
    get_edge(source,target)
      Returns the igraph edge going from source to target
    
    get_edge_attr(source,target,attr)
      Returns attribute of igraph edge going from source to target
    
    remove_vertices(_remove)
      Remove vertices from graph in list _remove
    
    assign_membership(membership,key)
      Assign vertices membership to some group. key defines the name
      of the membership. Vertices can be assigned to multiple memberships
      Membrship can be a list or dictionary
    
    assign_membership_list(membership,key)
      Use if membership is a list
    
    assign_membership_dict(membership,key)
      Use if membership is a dictionary
    
    group_vertices(key,combine_attrs='first')
      Group vertices beased on membership key.
    
    get_numpy_array(directed=False,vertex_order=None,edge_attr=None)
      Return adjacency matrix as a numpy array. 
    
    symmetrize_matrix()
      Will convert directed graph to undirected
    
    get_neighbors(vertex,mode="OUT",vattr="name")
      Return the neighbors of vertex
    
    map_vertex_names(vmap)
      Changes vertex names to names in dictionary vmap
    
    compute_strength(weight='weight')
      Computes normalized values of the edge weights

    reduce_to(G)
      Removes edges in graph not in graph G

    threshold_edge_greater_than(eattr,val)
      Removes edges with edge attribute less than or equal to val

    threshold_edge_less_than(eattr,val)
      Remove edges with edge attributes greater than or equal to val
      

    """
    def __init__(self,directed=False):
        """
        Parameters:
        -----------
        directed : bool (default False)
          If true, graph will be directed.
        
        """
        
        igraph.Graph.__init__(self,directed=directed)

    def __deepcopy__(self,memo):
        from copy import deepcopy
        return deepcopy(self)

    def get_edge(self,source,target):
        """
        Returns the igraph edge going from source to target
        
        Parameters:
        -----------
        source : str
           source vertex name
        target : str
           target vertex name
        
        Returns:
          igraph edge
        
        """
        
        try:
            i = self.vs.select(name=source)[0].index
            j = self.vs.select(name=target)[0].index
            return self.es.select(_source=i,_target=j)[0]
        except IndexError:
            return None

    def get_edge_attr(self,source,target,attr):
        """
        Returns attribute of igraph edge going from source to target
        
        Parameters:
        -----------
        source : str
          source vertex name
        target : str
          target vertex name
        attr : str
          attr key

        Returns:
        --------
        attribute value

        """

        i = self.vs.select(name=source)[0].index
        j = self.vs.select(name=target)[0].index
        w1 = self.es.select(_source=i,_target=j)[attr][0]
        #w2 = self.es.select(_source=j,_target=i)[attr][0]
        return w1
        
    def remove_vertices(self,_remove):
        """
        Remove vertices from graph in list _remove
        
        Parameters:
        ----------
        _remove : list
          List of vertex names to be removed

        """
        for n in _remove:
            for v in self.vs(name=n):
                v.delete()
                
    def assign_membership(self,membership,key='member'):
        """
        Assign vertices membership to some group. key defines the name
        of the membership. Vertices can be assigned to multiple memberships
        Membrship can be a list or dictionary   
        
        Parameters:
         membersip : list or dict
           Membership names for vertices
         key : str
           ID for the membership. Multiple memberships can be 
           assigned by using different membership keys
        
        """
    
        if isinstance(membership,list):
            self.assign_membership_list(membership,key=key)
        elif isinstance(membership,dict):
            self.assign_membership_dict(membership,key=key)

    def assign_membership_list(self,membership,key='member'):
        """
        Use if membership is a list
        
        Parameters:
         membersip : list
           Membership names for vertices. Order of list should 
           correspond to vertex indices.
         key : str
           ID for the membership. Multiple memberships can be 
           assigned by using different membership keys        
        
        """
        
        attr = [0]*self.vcount()
        for i in range(len(membership)):
            for j in membership[i]:
                attr[j] = i
        self.vs[key] = attr        

    def assign_membership_dict(self,membership,key='member'):
        """
        Use if membership is a dictionary
        
        Parameters:
         membersip : dict
           Membership names for vertices. Key is vertex name
           value is membership name.
         key : str
           ID for the membership. Multiple memberships can be 
           assigned by using different membership keys        
        
        """        

        attr = [0]*self.vcount()
        for v in self.vs:
            if v['name'] in membership:
                attr[v.index] = membership[v['name']]
            else:
                attr[v.index] = v[key]
        self.vs[key] = attr
        
    def group_vertices(self,key,combine_attrs='first'):
        """
        Group vertices based on membership key. 
        
        Parameters:
        -----------
        key : str
          Key ID for the vertex membership.
        combine_attrs : str
          Determines how to combine vertex attributes. See igraph API.
        
        """
        
        if isinstance(self.vs[key][0],str):
            vals = list(set(self.vs[key]))
            imap = dict([(vals[i],i) for i in range(len(vals))])
            membership = [imap[v[key]] for v in self.vs]
            self.contract_vertices(membership,combine_attrs=combine_attrs)
        else:
            self.contract_vertices(self.vs[key],combine_attrs=combine_attrs)
        self.simplify(loops=False,combine_edges=sum)
                
    def get_numpy_array(self,directed=False,vertex_order=None,edge_attr=None):
        """
        Return adjacency matrix as a numpy array.
        
        Paramters:
        ----------
        directed : bool
          If false, then the array will be symmetric
        vertex_order : list
          List of vertex names used to order rows and columns
        edge_attr : str
          Id of eddge attribute to use the matrix entries, e.g. the edge 'weights'
        
        Returns:
        --------
        Numpy array

        """
        
        if vertex_order:
            n = len(vertex_order)
            vhash = dict([(vertex_order[i],i) for i in range(n)])
            A = np.zeros((n,n))
            for e in self.es:
                w = 1
                if edge_attr: w = e[edge_attr]
                A[vhash[e.source]][vhash[e.target]] = w
                if not directed:
                    A[vhash[e.target]][vhash[e.source]] = w
        else:
            n = self.vcount()
            A = np.zeros((n,n))
            for e in self.es:
                w = 1
                if edge_attr: w = e[edge_attr]
                A[e.source][e.target] = w
                if not directed:
                   A[e.target][e.source] = w 
                
        return A

    def symmetrize_matrix(self):
        """
        Will convert directed graph to undirected
        """
        self.to_undirected(mode="collapse",combine_edges=sum)
        
    def get_neighbors(self,vertex,mode='OUT',vattr='name'):
        """
        Return the neighbors of vertex
        
        Parameters:
        -----------
        vertex : str
          Vertex name
        mode : str
          Edge direction for neighbors. 'IN' and 'OUT' choices.
          IN will look and indegree neighbors. OUT is the out degree neighbors.
        vattr : str
          Specifies which vertex attribute to return. 
        
        Return:
        -------
          Neighbor vertex attribute.

        """
        _neigh = self.neighbors(vertex, mode=mode)
        return self.vs[_neigh][vattr]
        
    def map_vertex_names(self,vmap):
        """
        Changes vertex names to names in dictionary vmap
        
        Parameters:
        ----------
        vmap : dict
         Dictionay used to map vertex names: 
         (key,val) = (old_vertex_name,new_vertex_name)
        
        Returns:
        --------
          Graph with new vertex name
        
        """
        G = Network(directed=self.is_directed())
        G.add_vertices([v['name'] for v in self.vs])
        for (i,e) in enumerate(self.es):
            G.add_edge(e.source,e.target)
            G.es[i]['weight'] = e['weight']
            G.es[i]['count'] = e['count']
        
        for v in G.vs:
            v['name'] = vmap[v['name']]
        
        return G

    def compute_strength(self,weight='weight'):
        """
        Computes normalized values of the edge weights

        Parameters:
        -----------
        weight : str
          Edge attribute id
        
        """
        v = range(self.vcount())
        self.vs['out-strength'] = self.strength(v,mode='out',
                                                loops=False,weights=weight)
        self.vs['in-strength'] = self.strength(v,mode='in',
                                               loops=False,weights=weight)
        for e in self.es:
            w = e['weight']
            sout = float(self.vs[e.source]['out-strength'])
            sin = float(self.vs[e.target]['in-strength'])
            e['out-strength'] = w / sout
            e['in-strength'] = w / sin

           
        
    def reduce_to(self,G):
        """
        Removes edges in graph not in graph G

        Parameters:
        -----------
        G : igraph
         Reference graph. Only edges in this graph will be maintained.
        
        """
        eid = []
        mode  = 0
        if not self.is_directed() and G.is_directed():
            mode = 1
        if mode == 1:
            for e in self.es:
                itoj = G.are_connected(e.source,e.target)
                jtoi = G.are_connected(e.target,e.source)
                if not itoj and not jtoi:
                    eid.append(e.index)
        else:
            for e in self.es:
                if not G.are_connected(e.source,e.target):
                    eid.append(e.index)
        self.delete_edges(eid)


    def threshold_edge_greater_than(self,eattr,val):
        """
        Removes edges with edge attribute less than or equal to val
        
        Parameters:
        -----------
        eattr : str
          Edge attribute name
        val : int or float
          Threshold value

        """
        es = [e for e in self.es if e[eattr] <= val]
        self.delete_edges(es)
        
    def threshold_edge_less_than(self,eattr,val):
        """
        Remove edges with edge attributes greater than or equal to val
        
        Parameters:
        -----------
        eattr : str
          Edge attribute name
        val : int or float
          Threshold value

        """        
        
        es = [e for e in self.es if e[eattr] >= val]
        self.delete_edges(es)
