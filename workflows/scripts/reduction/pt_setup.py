import networkx as nx 
from minigraphs.graph import spectral_radius

DICT_METRICS_FUNCS = {
    'density': nx.density,
    'assortativity_norm': lambda graph: (nx.degree_assortativity_coefficient(graph)+1)/2,
    'clustering': nx.average_clustering,
    'eig_1': lambda graph: spectral_radius(graph)
}