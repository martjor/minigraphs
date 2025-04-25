"""
Data Module
"""
import networkx as nx
from scipy.sparse import load_npz
from importlib.resources import path 

def load_graph(filename: str) -> nx.Graph:
    """Graph dataset.
    """
    filename = filename + '.npz'

    with path('minigraphs.data',filename) as file:
        return nx.from_scipy_sparse_array(load_npz(file))