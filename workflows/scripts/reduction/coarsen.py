from minigraphs import CoarseNET
from yaml import dump
from scripts.utils.io import load_graph, save_graph

# LOAD GRAPH
graph = load_graph(snakemake.input[0])

# COARSEN GRAPH
coarsener = CoarseNET(snakemake.params.alpha,
                      graph)
coarsener.coarsen()
graph_coarse = coarsener.G_coarse_

# SAVE GRAPH & METRICS
save_graph(snakemake.output[0], graph_coarse)

