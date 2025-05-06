import networkx as nx 
import numpy as np 
from networkx.linalg.algebraicconnectivity import algebraic_connectivity
from networkx.algorithms import edge_boundary
from scipy.sparse.linalg import eigsh

def edge_cut_size(G_full: nx.Graph, sub_nodes) -> int:
    """
    Number of edges with exactly one endpoint in the subgraph defined by `sub_nodes`
    and the other outside.

    Parameters
    ----------
    G_full : nx.Graph
        The original large graph (needed to see external edges).
    sub_nodes : Iterable
        The node set of your miniature (can be a set, list, or view).

    Returns
    -------
    int
        |∂S|  -- the size of the edge boundary of S.
    """
    # edge_boundary returns an *iterator* of edges → just count them
    return sum(1 for _ in edge_boundary(G_full, sub_nodes))

def graph_spectrum(graph: nx.Graph, k: int=1) -> float:
    """Calculates the spectral radius of a graph.
    """
    evals, _ = eigsh(nx.to_scipy_sparse_array(graph, dtype=np.float32), k=k, which='LM')
    
    return evals[0]

# Requires SciPy for sparse eigensolver; falls back to dense if small.
def laplacian_connectivity(G: nx.Graph, normalized: bool = False) -> float:
    """
    Compute the algebraic connectivity λ₂(L) of graph G.

    Parameters
    ----------
    G : nx.Graph
    normalized : bool
        If True, uses the normalized Laplacian.

    Returns
    -------
    float
        λ₂(L). Returns 0.0 for a graph with <2 nodes.
    """
    if G.number_of_nodes() < 2:
        return 0.0
    
    return float(algebraic_connectivity(G, normalized=normalized))