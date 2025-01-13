# Functions associated with manipulation of a graph

import networkx as nx 
import scipy as sp
import numpy as np

def spectral_radius(graph):
    '''Calculates the spectral radius of the adjacency matrix associated with a graph.
    '''
    adjacency = nx.adjacency_matrix(graph)._asfptype()
    evals,_ = sp.sparse.linalg.eigs(adjacency,k=1)
    return float(np.real(evals[0]))


