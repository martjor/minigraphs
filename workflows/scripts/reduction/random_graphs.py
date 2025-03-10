from scripts.utils.io import load_dict, save_graph
import networkx as nx 

# Load the dominant eigenvalue
eig = load_dict(snakemake.input[0])['eig_1']

# In an ER graph, the mean degree approximates the dominant eigenvalue.
# Estimate p from the number of edges and the eigenvalue.
p = eig / (snakemake.params.n_nodes-1)

# Construct graph
graph = nx.erdos_renyi_graph(snakemake.params.n_nodes, p)

# Save graph
save_graph(snakemake.output[0], graph)

