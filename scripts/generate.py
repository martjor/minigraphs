import networkx as nx 
from scipy.sparse import save_npz

generator = {
    "ER": nx.erdos_renyi_graph 
}

# Generate a complete graph
name = snakemake.params.params_gen['name']
n_vertices = snakemake.params.params_gen['n_vertices']
params = snakemake.params.params_gen['params'].values()

graph = generator[name](n_vertices,*params)

# Save the graph
save_npz(snakemake.output[0],nx.to_scipy_sparse_array(graph))