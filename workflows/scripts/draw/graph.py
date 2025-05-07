from numpy import load
import matplotlib
import matplotlib.pyplot as plt
from scripts.utils.io import load_graph, load_dict
import networkx as nx

matplotlib.use('Agg')
# Load graph and layout
graph = load_graph(snakemake.input[0])
layout = load(snakemake.input[1],allow_pickle=True)[()]

# Modify properties
properties = snakemake.params.properties
properties['edges']['width'] = 1.0
properties['nodes']['node_color'] = nx.clustering(graph).values()

# Draw graph with specified parameters
fig, ax = plt.subplots(figsize=(10,10), dpi=300)

nx.draw_networkx_edges(
    graph,
    layout,
    **properties['edges']
)

nx.draw_networkx_nodes(
    graph,
    layout,
    **properties['nodes']
)

ax.set_axis_off()

plt.savefig(snakemake.output[0], transparent=True, bbox_inches='tight')
