from utils import load_graph, save_dict
from networkx import write_gexf
from simulate import Contagion
import numpy as np

graph = load_graph(snakemake.input[0])
eq = np.load(snakemake.input[1])
n_trials = snakemake.params.n_trials
n_iterations = snakemake.params.n_iterations

# Calculate cdfs
cdf = np.cumsum(eq,axis=1)

nodes_red = np.searchsorted(cdf[0],np.random.rand(n_trials))
nodes_blue = np.searchsorted(cdf[1],np.random.rand(n_trials))

# Instantiate simulation object
count = {'red':0, 'blue':0}
for i in range(n_trials):
    simulation = Contagion(graph,
                        nodes_red[i],
                        nodes_blue[i],
                        n_iterations)
    
    simulation.run()

    count['red'] += simulation.count['red'] / n_trials
    count['blue'] += simulation.count['blue'] / n_trials

# Save data
save_dict(snakemake.output[0],count)
write_gexf(simulation.graph,snakemake.output[1])