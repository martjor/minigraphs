import numpy as np
from scripts.utils.io import load_dict, load_graph, save_graph

# Load graph and metrics
graph = load_graph(snakemake.input[0])
metrics = load_dict(snakemake.input[1])

# Calculate number of edges to remove
n_edges = graph.number_of_edges() - metrics['n_edges']

# Sparsify graph
edges = list(graph.edges)

# Remove edges at random
choice = np.random.choice(len(edges),size=n_edges,replace=False)
edges = [edges[idx] for idx in choice]
graph.remove_edges_from(edges)

# Save graph
save_graph(snakemake.output[0],graph)