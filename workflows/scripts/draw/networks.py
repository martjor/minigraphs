import numpy as np
import matplotlib.pyplot as plt
from scripts.utils.io import load_graph, load_dict
import networkx as nx

# Compute the number of miniatures
n_graphs = snakemake.params.n_graphs
n_miniatures = int(len(snakemake.input.miniatures) / (n_graphs * 3))

min_size = (5,5)
height = min_size[1] * n_miniatures
width = height + min_size[0] * n_graphs

# Create figure
fig = plt.figure(
    layout=None,
    figsize=(width, height),
    dpi=400
)

# Create custom grid
gs = fig.add_gridspec(
    nrows=n_miniatures,
    ncols=n_miniatures+n_graphs,
)

properties = {
    'nodes':{
        'node_size': 10,
        'node_color': '#18453B',
        'alpha': 0.05
    },
    'edges':{
        'width': .5,
        'edge_color': 'k',
        'alpha': 0.8
    }
}

def draw_graph(adj_file, metrics_file, layout_file, properties, ax):
    # Load files
    graph = load_graph(adj_file)
    metrics = {metric: value for metric, value in load_dict(metrics_file).items() if metric in (snakemake.params.targets + ['eig_1'])}
    layout = np.load(layout_file, allow_pickle=True)[()]
    
    # Draw graph
    nx.draw_networkx_nodes(
        graph, 
        layout,
        **properties['nodes']
    )
    
    # Draw edges
    nx.draw_networkx_edges(
        graph, 
        layout,
        **properties['edges']
    )
    
    text = '\n'.join([f"{metric}:{value:.4f}" for metric, value in metrics.items()])
    ax.text(
        0.0,
        1.0,
        text,
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes,
        color='blue'
    )

# Create axes for subraphs
axes = np.zeros((n_miniatures,n_miniatures+n_graphs),dtype=object)
for i in range(n_miniatures):
    for j in range(n_graphs):
        axes[i,j] = fig.add_subplot(gs[i, n_miniatures+j])
        
        # Obtain files
        idx = i * n_graphs * 3 + j * n_graphs
        files = snakemake.input.miniatures[idx:idx+3]
        
        draw_graph(
            *files,
            properties,
            axes[i,j]
        )
        
        n_nodes = files[0].split('/')[3].split('_')[1]
        axes[i,j].set_title(f"{n_nodes} nodes")
        
# Create axes for original graph
ax = fig.add_subplot(
    gs[:n_miniatures,:n_miniatures]
)
ax.set_title('Original')

properties['nodes']['alpha'] = 0.05
properties['edges']['alpha'] = 0.05
draw_graph(*snakemake.input.original, properties, ax)

fig.suptitle(f"{snakemake.wildcards.network.capitalize()} Network and Miniatures", fontsize=30)
plt.savefig(snakemake.output[0])