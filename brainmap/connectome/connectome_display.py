"""
connectome_display.py

For conveniently handling connectome data in the python command shell.

Inherits connectome_networkx.Connectome, and thus uses networkx data structures.

@author Christopher Brittin
@data 26 February 2019

"""
from tabulate import tabulate

from connectome.connectome_networkx import Connectome as _Connectome

class Connectome(_Connectome):
    """
    Class to handle connectome data in python command shell.

    Inherits connectome_networkx.Connectome 

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
    -----------
    ** see connectome_networkx.Connectome

    list_adjacency_neighbors(cell,attr="weigth")
        Lists the adjacency neighbors in order of attr
    
    list_chemical_neighbors(cell,attr="weigth")
        Lists the adjacency neighbors in order of attr

    list_electrical_neighbors(cell,attr="weigth")
        Lists the adjacency neighbors in order of attr
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
        _Connectome.__init__(self,db,neurons)

    def list_adjacency_neighbors(self,cell,attr='weight'):
        """
        Lists the adjacency neighbors in order of attr

        Parameters
        -----------
        cell : str
            Cell name 
        attr : str, default='weight'
            Edge attribute to sort neighbors
    
        """
        neigh = dict([(n,self.A[cell][n][attr]) for n in self.A.neighbors(cell)])
        sneigh = [(k, neigh[k]) for k in sorted(neigh, key=neigh.get, reverse=True)]
        print(tabulate(sneigh,headers=['Neighbor',attr]))

    def list_chemical_neighbors(self,cell,attr='weight'):
        """
        Lists the chemical neighbors in order of attr

        Parameters
        -----------
        cell : str
            Cell name 
        attr : str, default='weight'
            Edge attribute to sort neighbors
    
        """
        neigh = dict([(n,self.C[cell][n][attr]) for n in self.C.neighbors(cell)])
        sneigh = [(k, neigh[k]) for k in sorted(neigh, key=neigh.get, reverse=True)]
        print(tabulate(sneigh,headers=['Neighbor',attr]))

    def list_electrical_neighbors(self,cell,attr='weight'):
        """
        Lists the electrical neighbors in order of attr

        Parameters
        -----------
        cell : str
            Cell name 
        attr : str, default='weight'
            Edge attribute to sort neighbors
    
        """
        neigh = dict([(n,self.E[cell][n][attr]) for n in self.E.neighbors(cell)])
        sneigh = [(k, neigh[k]) for k in sorted(neigh, key=neigh.get, reverse=True)]
        print(tabulate(sneigh,headers=['Neighbor',attr]))
