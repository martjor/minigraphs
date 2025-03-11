import networkx as nx 
from scripts.utils.io import save_graph
import matplotlib.pyplot as plt

eig_1 = 6.9
def parameters_er(n_nodes):
    # Calculate optimal p
    p = eig_1 / n_nodes
    
    return nx.erdos_renyi_graph(n_nodes, p)

def parameters_ws(n_nodes):
    return nx.watts_strogatz_graph(n_nodes, round(eig_1), 0.005)

def parameters_ba(n_nodes):
    parameters = [2304, 576, 256]
    
    return nx.barabasi_albert_graph(n_nodes, parameters.index(n_nodes) + 1)

# Generate graphs for each generator
for file in snakemake.output:
    model, n_nodes = file.split('/')[-2].split('_')
    
    graph = eval(f"parameters_{model}({n_nodes})")
    
    save_graph(file,graph)
