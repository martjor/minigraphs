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
    
    Parameters
    ----------
    graph : nx.Graph
        Graph to sample from.
    n_nodes : int
        Number of nodes in the subgraph.
    n_swaps : int 
        Number of nodes to be exchanged from the boundary at each iteration.
    seed : int 
        State of the random number generator.
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

class RandomizedSet:
    def __init__(self, random_state: random.Random):
        self.val_to_index = {}
        self.values = []
        self.random = random_state

    def insert(self, val: int) -> bool:
        if val in self.val_to_index:
            return False
        self.val_to_index[val] = len(self.values)
        self.values.append(val)
        return True

    def remove(self, val: int) -> bool:
        if val not in self.val_to_index:
            return False
        idx = self.val_to_index[val]
        last_val = self.values[-1]
        self.values[idx] = last_val
        self.val_to_index[last_val] = idx
        self.values.pop()
        del self.val_to_index[val]
        return True

    def getRandom(self) -> int:
        return self.random.choice(self.values)
    
    def __bool__(self):
        return bool(self.values)

class CoarseningChain(Chain):
    """A graph chain implementing node coarsening.
    """

    def __init__(self, graph: nx.Graph, n_groups: int, seed: int=None):
        super().__init__(seed=seed)
        self.graph = graph 
        self.n_groups = n_groups
        self.groups = range(n_groups)

        # Generate coarsened graph.
        # Nodes hold references to sets containing the nodes encapsulated by the super nodes
        # Edges contain scalars reflecting the number of edges a super edge represents
        self.graph_coarse: nx.Graph = nx.complete_graph(self.n_groups)
        nx.set_node_attributes(self.graph_coarse, values = {node: {'nodes': RandomizedSet(self.random)} for node in self.graph_coarse.nodes()})
        nx.set_edge_attributes(self.graph_coarse, name='count', values=0)

        # Randomly allocate nodes to super nodes
        for node in self.graph.nodes():
            # Select a group at random
            group = self.random.choice(self.groups)

            # Store references
            self.graph.nodes[node]['group'] = group 
            self.graph_coarse.nodes[group]['nodes'].insert(node)

        # Identify edges of coarsened graph
        for edge in self.graph.edges():
            # Identify corresponding super edge
            super_edge = tuple((self.graph.nodes[node]['group'] for node in edge))

            if super_edge[0] != super_edge[1]:
                self.graph_coarse.edges[super_edge]['count'] += 1

    @property 
    def state(self) -> nx.Graph:
        """Returns the coarsened graph.
        """
        graph = nx.empty_graph(self.n_groups)
        graph.add_edges_from(((u,v) for (u,v,data) in self.graph_coarse.edges(data=True) if data['count'] > 0))
        
        return graph
    
    def propose(self) -> nx.Graph:
        """Proposes a new graph by randomly moving nodes across super nodes.

        Randomly selects a super node to draw a node and moves it to a random different super node,
        updating the edge counts along the way.
        """
        # Select random group
        old_group = self.random.choice(self.groups)
        nodes: RandomizedSet = self.graph_coarse.nodes[old_group]['nodes']

        if nodes: # Check for nodes in group
            # Select random node
            node = nodes.getRandom()
            neighbors = list(self.graph.neighbors(node))

            if neighbors: # Check for neighbors
                # Select random neighbor
                neighbor = self.random.choice(neighbors)
                new_group = self.graph.nodes[neighbor]['group']

                if old_group != new_group:
                    # Update data.
                    self.change = (node, old_group)
                    self._move_node(node, new_group)

    def _move_node(self, node: int, new_group: int):
        """Moves a node to a new group.
        """
        old_group = self.graph.nodes[node]['group']

        # Update node data
        self.graph.nodes[node]['group'] = new_group
        self.graph_coarse.nodes[old_group]['nodes'].remove(node)
        self.graph_coarse.nodes[new_group]['nodes'].insert(node)

        # Update edge counts
        for neighbor in list(self.graph.neighbors(node)):
            neighbor_group = self.graph.nodes[neighbor]['group']

            # Update edge counts
            if old_group != neighbor_group:
                old_edge = (old_group, neighbor_group)
                self.graph_coarse.edges[old_edge]['count'] -= 1 

            if new_group != neighbor_group:
                new_edge = (new_group, neighbor_group)
                self.graph_coarse.edges[new_edge]['count'] += 1

    def reject(self):
        """Reject the proposed change and reverse the proposed merge.
        """
        # Update node data
        node, old_group = self.change
        self._move_node(node, old_group)
