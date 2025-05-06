from typing import Callable, Iterable, Union, Optional, Deque
import math
import random
from tqdm import tqdm
import networkx as nx
import numpy as np
from abc import ABC, abstractmethod
from collections import deque

class Chain(ABC):
    """Base class for graph‑valued Markov chains.

    Parameters
    ----------
    graph : nx.Graph
        The reference graph from which subgraphs will be proposed.
    current_graph : Optional[nx.Graph]
        Initial state of the chain. If ``None``, use ``original_graph``.
    rng : Optional[random.Random]
        RNG instance so that experiments can be reproduced.
    """

    def __init__(
        self,
        graph: nx.Graph,
        seed: Optional[int] = None,
    ) -> None:
        self.graph = graph
        self.random = random.Random(seed)
        
    @abstractmethod
    def _propose(self) -> nx.Graph:
        """Proposes a new graph from the current state
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

class SubGraphChain(Chain):
    def __init__(
            self, 
            graph,
            n_nodes,
            n_swaps=1,
            seed = None
        ) -> None:
        super().__init__(graph, seed)
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

class SimulatedAnnealing:
    """Generic simulated‑annealing / MCMC driver on top of a ``Chain``.

    Parameters
    ----------
    chain : Chain
        Instance that supplies proposals via ``chain.propose()``.
    energy : Callable[[nx.Graph], float]
        Function that returns *energy* (lower is better) of a graph.
    schedule : Union[float, Callable[[int], float]]
        Temperature schedule. A constant (float) means fixed temperature.
        Otherwise provide a function ``T(step)``.
    burn_in : int, default 0
        Number of initial iterations to discard.
    thinning : int, default 1
        Keep one sample every ``thinning`` accepted states.
    max_steps : int, default 10_000
        Total number of proposal steps.
    rng : Optional[random.Random]
        Random generator (for reproducibility).
    """

    def __init__(
        self,
        chain: Chain,
        energy: Callable[[nx.Graph], float],
        schedule: Union[float, Callable[[int], float]],
        *,
        max_steps: int = 1_000,
        seed: Optional[int] = None,
    ) -> None:
        self.chain = chain
        self.energy_fn = energy
        self.schedule = schedule if callable(schedule) else lambda _: float(schedule)
        self.max_steps = max_steps
        self.random = random.Random(seed)

        # Book‑keeping
        self.samples: Deque[nx.Graph] = deque()
        self.energies: Deque[float] = deque()
        self.acceptance_count = 0
        self.total_proposals = 0

    def run(self) -> None:
        """Run the annealing / MCMC loop."""
        graph_curr = self.chain.state
        energy_curr = self.energy_fn(graph_curr)

        for step in tqdm(range(self.max_steps)):
            beta = self.schedule(step)
            graph_new = self.chain._propose()
            energy_new = self.energy_fn(graph_new)

            # Metropolis acceptance probability
            dE = energy_new - energy_curr
            accept = dE < 0 or self.random.random() < math.exp(-beta * dE)

            self.total_proposals += 1
            if accept:
                self.acceptance_count += 1
                graph_curr, energy_curr = graph_new, energy_new
                self.chain.state = graph_curr

            self.samples.append(graph_curr)
            self.energies.append(energy_curr)

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------
    @property
    def acceptance_rate(self) -> float:
        return self.acceptance_count / max(1, self.total_proposals)

    def summary(self) -> str:
        return (
            f"SimulatedAnnealing: {len(self.samples)} samples | "
            f"acc_rate={self.acceptance_rate:.3f} | "
            f"energy_mean={np.mean(self.energies) if self.energies else float('nan'):.3f}"
        )

    # Convenience iterator
    def __iter__(self) -> Iterable[nx.Graph]:
        return iter(self.samples)
