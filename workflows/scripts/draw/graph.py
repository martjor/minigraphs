from numpy import load
import matplotlib.pyplot as plt
from scripts.utils.io import load_graph, load_dict
import networkx as nx
from module import draw_graph

# Load graph and layout
graph = load_graph(snakemake.input[0])
layout = load(snakemake.input[1],allow_pickle=True)[()]

# Draw graph with specified parameters
fig, ax = plt.subplots(figsize=(10,10), dpi=300)

draw_graph(
    graph,
    layout,
    snakemake.params.properties
)

plt.savefig(snakemake.output[0])
