import random
import networkx as nx
from abc import ABC, abstractmethod

class Chain(ABC):
    """Base class for graph‑valued Markov chains.

    Parameters
    ----------
    seed : int
        Random state of the chain.
    """

    def __init__(
        self,
        seed: int=None,
    ) -> None:
        self.random = random.Random(seed)
        
    @abstractmethod
    def _propose(self) -> nx.Graph:
        """Proposes a new graph from the current state.
        """
    
    def __iter__(self):
        return self 
    
    def __next__(self):
        """Propse new graph and update internal state.
        """
        graph = self._propose()
        self.graph_current = graph 
        return graph

    @property
    def state(self) -> nx.Graph:
        return self.graph_current
    
    @state.setter
    def state(self, graph: nx.Graph):
        self.graph_current = graph

class SubgraphUniform(Chain):
    """Subgraphs are proposed by sampling uniformly from the original pool of nodes.
    """
    def __init__(self, graph, n_nodes, seed=None):
        super().__init__(seed)
        self.graph = graph
        self.n_nodes = n_nodes

        # Random subgraph
        self.graph_current = nx.subgraph(graph, self.random.sample(list(graph.nodes), k=self.n_nodes))

    def _propose(self):
        """Proposes a uniform random subgraph
        """
        return nx.subgraph(self.graph, self.random.sample(list(self.graph.nodes), k=self.n_nodes))

class SubgraphBoundary(Chain):
    """Proposes subgraphs by replacing `n_swaps` nodes for nodes on the boundary with the subgraph. 
    """
    def __init__(
            self, 
            graph,
            n_nodes,
            n_swaps=1,
            seed = None
        ) -> None:
        super().__init__(seed)
        self.graph = graph
        self.n_nodes = n_nodes
        self.n_swaps = n_swaps

        # Generate subgraph in the largest connected component
        nodes = list(max(nx.connected_components(graph), key=len))

        if len(nodes) <= n_nodes:
            raise ValueError("Largest connected component is not large enough.")
        else:
            self.graph_current = nx.subgraph(graph, self.random.sample(nodes, k=self.n_nodes))

    def _propose(self):
        """Proposes a new graph by swapping nodes from the boundary.
        """
        nodes = list(self.graph_current.nodes)
        nodes_boundary = list(nx.node_boundary(self.graph, nodes))

        nodes_new = (
            self.random.sample(nodes, k=(self.n_nodes - self.n_swaps)) + 
            self.random.sample(nodes_boundary, k=self.n_swaps)
        )

        return nx.subgraph(self.graph, nodes_new)

