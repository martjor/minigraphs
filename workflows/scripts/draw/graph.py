import matplotlib.pyplot as plt 
import networkx as nx
from scripts.utils.io import load_graph
from numpy import load

# Load the graph and Layout
graph = load_graph(snakemake.input[0])
pos = load(snakemake.input[1],allow_pickle=True)[()]

# Draw Graph
fig, ax = plt.subplots(figsize=(5,5),dpi=600)
color = list(nx.clustering(graph).values())

ax.set_title(f"|V|={graph.number_of_nodes()}, den={nx.density(graph):.2e}")
ax.set_yticks([])
ax.set_xticks([])
ax.set_facecolor('black')

nx.draw_networkx_nodes(graph,pos,
                       node_size=0.1,
                       node_color=color,
                       vmin=0.0,
                       vmax=1.0,
                       alpha=0.7,
                       cmap='cool')
nx.draw_networkx_edges(graph,pos,
                       edge_color="white",
                       alpha=0.05,
                       width=0.5)

plt.savefig(snakemake.output[0],bbox_inches='tight')