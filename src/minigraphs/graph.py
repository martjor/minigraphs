# Functions associated with manipulation of a graph

import networkx as nx 
import scipy as sp
import numpy as np
from functools import lru_cache
from scipy.stats import moment

def spectral_radius(graph):
    '''Calculates the spectral radius of the adjacency matrix associated with a graph.
    '''
    adjacency = nx.adjacency_matrix(graph)._asfptype()
    evals,_ = sp.sparse.linalg.eigs(adjacency,k=1)
    return float(np.real(evals[0]))

@lru_cache(maxsize=5)
def degree_sequence(graph):
    '''Computes the node degree sequence of a graph
    '''
    sequence = sorted((d for n, d in graph.degree()))
    return sequence 

@lru_cache(maxsize=5)
def degree_distribution_moment(graph,order: int = 1,center: float = 0.0):
    sequence = degree_sequence(graph)
    '''Computes the specified moment from the degree distribution of the graph
    '''
    
    return float(moment(sequence,order,center=center))


