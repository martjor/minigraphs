import networkx as nx 
from scripts.utils.io import save_graph

eig_1 = 14.1421356237
def parameters_er(n_nodes):
    # Calculate optimal p
    p = eig_1 / n_nodes
    
    return nx.erdos_renyi_graph(n_nodes, p)

def parameters_ws(n_nodes):
    return nx.watts_strogatz_graph(n_nodes, round(eig_1), 0.005)

def parameters_ba(n_nodes):    
    return nx.barabasi_albert_graph(n_nodes, snakemake.params.sizes.index(n_nodes) + 2)


graph = eval(f"parameters_{snakemake.wildcards.model}({snakemake.wildcards.n_nodes})")

save_graph(snakemake.output[0],graph)
