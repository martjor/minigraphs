import networkx as nx 
import matplotlib.pyplot as plt

def draw_graph(
        graph,
        layout,
        properties,
    ) -> None:
    '''Draws a graph using the specified layout.
    
    Draws the graph and allows for granular control of the appearance  of nodes and edges
    '''
    
    # Draw Nodes
    nx.draw_networkx_nodes(
        graph, 
        layout,
        **properties['nodes'],
    )

    # Draw Edges
    nx.draw_networkx_edges(
        graph, 
        layout,
        **properties['edges'],
    )
    
    