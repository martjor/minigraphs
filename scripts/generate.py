import networkx as nx 
from scipy.sparse import save_npz

# Generate an ER graph with 10 vertices
graph = nx.erdos_renyi_graph(10,snakemake.params.p)

# Save the graph
save_npz(snakemake.output[0],nx.to_scipy_sparse_array(graph))