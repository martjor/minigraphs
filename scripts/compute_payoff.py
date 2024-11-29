from utils import load_graph, save_dict, save_graph
from simulate import Contagion

graph = load_graph(snakemake.input[0])
node_r = int(snakemake.wildcards.i)
node_b = int(snakemake.wildcards.j)
n_iterations = snakemake.params.n_iterations

# Instantiate simulation object
simulation = Contagion(graph,
                       node_r,
                       node_b,
                       n_iterations)

# Run simulation
simulation.run()

# Save data
save_dict(snakemake.output[0],simulation.count)
save_graph(snakemake.output[1],simulation.graph)

