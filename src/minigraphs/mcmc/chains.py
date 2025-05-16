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
    def propose(self):
        """Proposes a new graph and updates its internal state.
        """
        pass 

    @abstractmethod
    def reject(self):
        """Rejects the graph proposal and retrieves the previous internal state.
        """
        pass

    @property
    @abstractmethod
    def state(self) -> nx.Graph:
        """Retrieves the current state of the chain.
        """
        
    def __iter__(self):
        return self 
    
    def __next__(self):
        """Propse new graph and update internal state.
        """
        state = self.state 
        self.propose()
        return state

class SubgraphUniform(Chain):
    """Subgraphs are proposed by sampling uniformly from the original pool of nodes.
    """
    def __init__(self, graph, n_nodes, seed=None):
        super().__init__(seed)
        self.graph = graph
        self.n_nodes = n_nodes

        # Random subgraph
        self.subgraph = nx.subgraph(graph, self.random.sample(list(graph.nodes), k=self.n_nodes))

    def propose(self):
        """Proposes a uniform random subgraph
        """
        self.old_nodes = [node for node in self.subgraph.nodes()]
        self.subgraph = nx.subgraph(self.graph, self.random.sample(list(self.graph.nodes), k=self.n_nodes))
    
    def reject(self):
        self.subgraph = self.subgraph(self.graph, self.old_nodes)
    
    @property 
    def state(self):
        """The currently sampled subgraph
        """
        return self.subgraph
    


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
            self.subgraph = nx.subgraph(graph, self.random.sample(nodes, k=self.n_nodes))

    def propose(self):
        """Proposes a new graph by swapping nodes from the boundary.
        """
        self.old_nodes = list(self.subgraph.nodes())
        nodes_boundary = list(nx.node_boundary(self.graph, self.old_nodes))

        nodes_new = (
            self.random.sample(self.old_nodes, k=(self.n_nodes - self.n_swaps)) + 
            self.random.sample(nodes_boundary, k=self.n_swaps)
        )

        self.subgraph = nx.subgraph(self.graph, nodes_new)
    
    def reject(self):
        self.subgraph = nx.subgraph(self.graph, self.old_nodes)
    
    @property 
    def state(self):
        """The currently sampled subgraph
        """
        return self.subgraph

