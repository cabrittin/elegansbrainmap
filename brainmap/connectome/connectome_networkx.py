"""
connectome_networkx.py

Connectome data structures. Uses Networkx 

See https://networkx.github.io/

Author: Christopher Brittin

"""

import networkx as nx
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
    -------
    update_neurons(_neurons)
      Update the neuron list
    
    remove_cells(vertices)
      Removes vertices from graphs C, E and A.

    remove_vertices(G,remove)
      Remove vertices from graph G
    
    group_cells(groups,**kwargs)
      Group cells based on dictionary groups; 
      (key,value) = (cell_name,group_name)

    group_vertices(G,GROUPS)
      Group vertrices in graph G based on dictionary GROUPS

    load_chemical(synapses,add_poly=False)
      Create chemical connectivity graph by loading edges from synapses. 
      If add_poly, then the number of polyad and monad synapses is tracked.

    load_electrical(synapses)
      Create gap junction connectivity graph by loading edges from synapses. 
      If add_poly, then the number of polyad and monad synapses is tracked.    

    load_edges(G,vertices,edges,add_poly=False)
      Load edges between vertices into graph G. If add_poly, then the
      number of polyad and monad synapses is tracked.    

    load_adjacency(adjacency,directed=False)
      Load adjacency graph from _adjacency edges. If directed, adjacency graph
      will be made directed. 
    
    reduce_to_adjacency()
      Reduce chemical and gap junction connectivity graphs to nodes and
      edges found in the adjacency graph

    _reduce_to_adjacency(A,H)
      Remove edges in graph H that are not in graph A

    combine_chem_and_elec()
      Combine the chemical and gap junction connectivity graphs

    remove_self_loop(self)
      Removes loops from graphs self.C and self.E

    _remove_self_loops(G)
      Removes loops from graph G

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
       
    def update_neurons(self,_neurons):
        """
        Update the neuron/cell/node list

        Parameters:
        -----------
        _neurons : list 
           list of neurons/cell/node/names

        """
        self.neurons = _neurons
        self.size = len(_neurons)

    def remove_cells(self,vertices):
        """"
        Remove vertices from graphs C, E and A
        
        Parameters
        ----------
        vertices : list
          List of vertex names

        """        
        if self.C:self.remove_vertices(self.C,vertices)
        if self.E:self.remove_vertices(self.E,vertices)
        if self.A:self.remove_vertices(self.A,vertices)

    def remove_vertices(self,G,remove):
        """
        Remove vertices in list remove from graph G
        
        Parameters
        ----------
        G : Networkx
         Graph from which vertices will be removed
        remove : list
         List of vertex names to be removed.
        
        """
        remove = set(remove) & set(G.nodes())
        for n in remove:
            G.remove_node(n)

    def group_cells(self,groups,**kwargs):
        """
        Group vertices based on dictionary groups. The grouping identified
        by key (default 'group'). So multiple groups can be assigned to 
        same graphs.                
       
        Parameters
        ----------
        groups : dict
          Dictionary of vertex memberships. 
          (key,value) = (cell_name,group_name)
        
        **kwargs : dummy variable 
          Not used in for Networkx but kept to remain consistent
          with igraph implementation
        """
        if self.C: self.C = self.group_vertices(self.C,groups)
        if self.E: self.E = self.group_vertices(self.E,groups)
        if self.A: self.A = self.group_vertices(self.A,groups)

    def group_vertices(self,G,GROUPS):
        """
        Group vertices based on dictionary GROUPS in graph G

        Parameters
        ----------
        G : Networkx
          Graph object
        GROUPS : dict
          Group assignments. (key,value) = (cell_name,group_name)
        
        """
        if nx.is_directed(G):
            H = nx.DiGraph()
        else:
            H = nx.Graph()

        for e in G.edges():
            attr = G[e[0]][e[1]]
            if e[0] in GROUPS:
                n1 = GROUPS[e[0]]
            else:
                n1 = e[0]
            if e[1] in GROUPS:
                n2 = GROUPS[e[1]]
            else:
                n2 = e[1]
            if (n1,n2) in H.edges():
                #print n1,n2,H[n1][n2]
                H[n1][n2]['weight'] += attr['weight']
                H[n1][n2]['count'] += attr['count']
            else:
                H.add_edge(n1,n2,weight = attr['weight'],
                           count=attr['count'])
        return H

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
        self.C = nx.DiGraph()
        self.load_edges(self.C,self.neurons,synapses,add_poly=add_poly)

    def load_electrical(self,synapses):
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
        
        self.E = nx.Graph()
        self.load_edges(self.E,self.neurons,synapses)
                
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
        
        for e in edges:
            pre = ioaux.format.rm_brack(e[0])
            if pre not in vertices: continue
            #i_pre = self.neurons[pre]
            _post = list(set(map(ioaux.format.rm_brack,e[1].split(','))))
            if len(_post) == 1:
                poly = 'S'
            else:
                poly = 'Sp'
            if self.db == 'N2U' and e[4] in ['VC','DC']:
                w = 2*int(e[2]) - 1
            else:
                w = int(e[2])
                
            for post in _post:
                if post not in vertices: continue 
                #i_post = self.neurons[post]
                if not G.has_edge(pre,post):
                    if add_poly:
                        G.add_edge(pre,post,weight=0.0,count=0,S=0,Sp=0)
                    else:
                        G.add_edge(pre,post,weight=0.0,count=0)
                G[pre][post]['weight'] += w
                G[pre][post]['count'] += 1
                if add_poly:
                    G[pre][post][poly] += 1

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
        
        if directed:
            self.A = nx.DiGraph()
        else:
            self.A = nx.Graph()
        for (i,j,weight,imgNum) in adjacency:
            weight = int(weight)
            count = 1
            if self.db == 'N2U' and 'VC' in imgNum:
                weight *=2
                count = 2
            if not self.A.has_edge(i,j):
                self.A.add_edge(i,j,weight=0,count=0)
            self.A[i][j]['weight'] += weight
            self.A[i][j]['count'] += count
            
    def reduce_to_adjacency(self):
        """
        Reduce chemical and gap junction connectivity graphs to nodes and
        edges found in the adjacency graph  
        """
        
        self.C = self._reduce_to_adjacency(self.A,self.C)
        self.E = self._reduce_to_adjacency(self.A,self.E)

    @staticmethod
    def _reduce_to_adjacency(A,H):
        """
        Eliminates nodes and edges in H not found in A
        
        Parameters:
        -----------
        A : Networkx
         Reference graph. Nodes and edges not found in A
         will be removed.
        G : Networkx
         Graph to be changed. Nodes and edges not found in A
         will be removed from H.
        
        """
        G = H.copy()
        GnotAnodes = list(set(G.nodes()) - set(A.nodes()))
        G.remove_nodes_from(GnotAnodes)
        edges = [e for e in G.edges()]
        for (a,b) in edges:
            if not A.has_edge(a,b):
                G.remove_edge(a,b)
        return G
            
                
    def combine_chem_and_elec(self):
        """
        Combine the chemical and gap junction connectivity graphs
        Combined graph stored in attribute self.D

        """
        self.D = self.C.copy()
        for (a,b) in self.E.edges():
            w = 0.5*self.E[a][b]['weight']
            if not self.D.has_edge(a,b): self.D.add_edge(a,b,weight=0)
            if not self.D.has_edge(b,a): self.D.add_edge(b,a,weight=0)
            self.D[a][b]['weight'] += w
            self.D[b][a]['weight'] += w

    def remove_self_loops(self):
        """
        Removes loops from graphs self.C and self.E
        """
        self.C = self._remove_self_loops(self.C)
        self.E = self._remove_self_loops(self.E)

    @staticmethod
    def _remove_self_loops(G):
        """
        Removes loops from graph G
        
        Parameters
        ----------
        G : Networkx
          Graph
        """
        if nx.is_directed(G):
            H = nx.DiGraph(G)
        else:
            H = nx.Graph(G)
        for (n1,n2) in G.edges():
            if n1 == n2: H.remove_edge(n1,n2)
        return H

        
    def split_left_right(self,left,right):
        """
        Split out left and right graphs
        
        Parameters:
        -----------
        left: list
            List of left cells
        right: list
            List of right cells
        """
        if self.A:
            self.Al = self.split_graph(self.A,left)
            self.Ar = self.split_graph(self.A,right)
        if self.C:
            self.Cl = self.split_graph(self.C,left)
            self.Cr = self.split_graph(self.C,right)
        if self.E:
            self.El = self.split_graph(self.E,left)
            self.Er = self.split_graph(self.E,right)


    @staticmethod
    def split_graph(G,nodes):
        """
        Splits out the graphs with all edges connected to given nodes

        Parameters:
        -----------
        G : networx Graph
            Graph
        nodes : list
          List of nodes to keep
        """
        if nx.is_directed(G):
            H = nx.DiGraph()
        else:
            H = nx.Graph()
        for n in nodes:
            if not G.has_node(n): continue
            for m in G.neighbors(n):
                H.add_edge(n,m,weight=G[n][m]['weight'])
            if G.is_directed():
                for m in G.predecessors(n):
                    H.add_edge(m,n,weight=G[m][n]['weight'])

        return H


    def map_right_graphs(self,nmap):
        """
        Maps the nodes in the right graphs to the left nodes

        Parameters:
        -----------
        nmap: dict
            dictionary to map right nodes to left nodes

        """
        if hasattr(self,'Ar'): self.Ar = self.map_graph_nodes(self.Ar,nmap)
        if hasattr(self,'Cr'): self.Cr = self.map_graph_nodes(self.Cr,nmap)
        if hasattr(self,'Er'): self.Er = self.map_graph_nodes(self.Er,nmap)

    @staticmethod
    def map_graph_nodes(G,lrmap):
        """
        Maps left and right nodes of graph

        G : networkx graph
            Graph
        lrmap : dict
            Dictionary that maps left to right nodes and vice versa
        """
        if nx.is_directed(G):
            H = nx.DiGraph()
        else:
            H = nx.Graph()

        for (a,b) in G.edges():
            H.add_edge(lrmap[a],lrmap[b],weight=G[a][b]['weight'])
        return H

