import matplotlib.pyplot as plt
import networkx as nx 
from typing import Optional, List, Union, Dict, Any, Tuple
from matplotlib.axes import Axes
from functools import wraps

__all__=[
    "draw_subgraph"
]

Graph = Union[nx.Graph, nx.DiGraph]


def get_context(func):
    """Decorator that wraps visualization functions to allow for automatic context retrieval.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Get current axes if non is specified
        ax = kwargs.get('ax')

        if ax is None:
            cf = plt.gcf()
        else:
            cf = ax.get_figure()

        cf.set_facecolor("w")
        if ax is None:
            if cf.axes:
                ax = cf.gca()
            else:
                ax = cf.add_axes((0, 0, 1, 1))

        kwargs['ax'] = ax 

        return func(*args, **kwargs)
    
    return wrapper

@get_context
def draw_subgraph(
    graph: Graph, 
    subgraph: Union[Graph, List[int]], 
    pos: Dict[int, Tuple[float,float]],
    * , 
    graph_kwargs: Optional[Dict[str, Any]]=None, 
    subgraph_kwargs: Optional[Dict[str, Any]]=None, 
    ax: Optional[Axes]=None
):
    """Draws a graph and highlights the nodes corresponding to the specified subgraph.

    Parameters
    ----------
    graph : Graph
        The original graph
    subgraph : Union[Graph, List[int]]
        The subgraph or set of nodes the finde the subgraph.
    pos : Dict[int, Tuple[float, float]]
        A dictionary with the node positions
    """
    # Retrieve nodelist
    nodelist = list(subgraph.nodes()) if isinstance(subgraph, Graph) else subgraph
    graph_kwargs = graph_kwargs if graph_kwargs is not None else {}
    subgraph_kwargs = subgraph_kwargs if subgraph_kwargs is not None else {}

    first_color = graph_kwargs.get('node_color')
    second_color = subgraph_kwargs.get('node_color')

    graph_kwargs['node_color'] = first_color if first_color is not None else 'tab:blue'
    subgraph_kwargs['node_color'] = second_color if second_color is not None else 'tab:orange'

    nx.draw_networkx_nodes(graph, pos, ax=ax, **graph_kwargs)
    nx.draw_networkx_nodes(graph, pos, nodelist=nodelist, ax=ax, **subgraph_kwargs)
    nx.draw_networkx_edges(graph, pos, ax=ax)

    return ax
    
