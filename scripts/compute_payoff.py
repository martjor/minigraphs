from utils import load_graph, save_dict
from networkx import write_gexf
from simulate import Contagion

graph = load_graph(snakemake.input[0])
node_r = int(snakemake.wildcards.i)
node_b = int(snakemake.wildcards.j)
n_iterations = snakemake.params.n_iterations
n_trials = snakemake.params.n_trials

# Instantiate simulation object
simulation = Contagion(graph,
                       node_r,
                       node_b,
                       n_iterations)

# Run simulation
count = {'red':0, 'blue':0}
for i in range(n_trials):
    simulation.run()

    count['red'] += simulation.count['red'] / n_trials
    count['blue'] += simulation.count['blue'] / n_trials

# Save data
save_dict(snakemake.output[0],count)
write_gexf(simulation.graph,snakemake.output[1])

