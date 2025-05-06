import networkx as nx 
import yaml
from numpy import save
from scripts.utils.io import load_graph

# Load graph
graph = load_graph(snakemake.input[0])

# Calculate Layout
pos = nx.forceatlas2_layout(
    graph,
    max_iter=snakemake.params.max_iter,
    dissuade_hubs=True,
    scaling_ratio=1.0,
)

# Save layout
save(snakemake.output[0],pos)

