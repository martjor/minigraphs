import networkx as nx 
import yaml
from numpy import save
from utils import load_graph

# Load graph
graph = load_graph(snakemake.input[0])

# Calculate Layout
pos = nx.forceatlas2_layout(graph)

# Save layout
save(snakemake.output[0],pos)

