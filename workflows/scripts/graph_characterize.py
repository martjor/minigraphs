from minigraphs.miniaturize import NX_ASSORTATIVITY, NX_CLUSTERING, NX_DENSITY
import networkx as nx 
from yaml import dump
from scipy.sparse import load_npz
from scipy.sparse.linalg import eigs
import numpy as np

def graph_components(G):
    # Get connected components of the graph
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    
    # Generate induced graphs
    graphs = [G.subgraph(c).copy() for c in components]

    return graphs

def degree_distribution(graph):
    '''Computes the node degree distribution of a graph
    '''
    hist = nx.degree_histogram(graph)
    dist = np.zeros((len(hist),2))
    dist[:,0] = np.arange(len(hist))
    dist[:,1] = hist
    
    # Normalize distribution
    dist[:,1] /= dist[:,1].sum()
    
    return dist 

def moment(dist,p):
    '''Calculates the pth momenth of a distribution
    '''
    moment = ((dist[:,0] ** p) * dist[:,1]).sum()
    return moment

def moment_central(dist,mean,p):
    '''Calculates the pth central moment of a distribution
    '''
    moment = np.power(np.abs(dist[:,1] * (dist[:,0] - mean)),p).sum()
    return moment 

# PREAMBLE 
adjacency_file = snakemake.input[0]
metrics_file = snakemake.output[0]

functions = {
    'n_nodes': lambda G: G.number_of_nodes(),
    'n_edges': lambda G: G.number_of_edges(),
    'density': NX_DENSITY,
    'clustering': NX_CLUSTERING,
    'assortativity': NX_ASSORTATIVITY,
}

# LOAD ADJACENCY MATRIX AND GENERATE GRAPH
adjacency = load_npz(adjacency_file)
adjacency = adjacency._asfptype()
graph = nx.from_scipy_sparse_array(adjacency)

# CALCULATE METRICS
components = graph_components(graph)

# Calculate eigenvalues
evals,_ = eigs(adjacency,2)

# Calculate degree distribution
dist = degree_distribution(graph)

# Evaluate metrics
metrics = {metric: func(graph) for metric, func in functions.items()}
metrics['n_components'] = len(components)
metrics['connectivity'] = components[0].number_of_nodes() / graph.number_of_nodes()
metrics['eig_1'] = float(np.real(evals[0]))
metrics['eig_2'] = float(np.real(evals[1]))
metrics['dd_m1'] = moment(dist,1)
metrics['dd_m2'] = moment(dist,2)

# WRITE METRICS
with open(metrics_file,'w') as file:
    dump(metrics,file,default_flow_style=False)
    