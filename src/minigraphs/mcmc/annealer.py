from typing import Union, Callable, Iterable, Optional, Deque
from collections import deque
from tqdm import tqdm 
from .chains import Chain 
from math import exp
import random 
import numpy as np 
from networkx import Graph

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
    max_steps : int, default 10_000
        Total number of proposal steps.
    seed : Optional[int]
        Random state for the ennealer.

    Attributes
    ----------
    energies_ : Deque[float]
        History of energies at every time step.
    best_graph_ : Tuple[Graph, float]
        Graph with the lowest energy throghout the process, along with it's corresponding energy.
    """
    def __init__(
        self,
        chain: Chain,
        energy: Callable[[Graph], float],
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

    def run(self) -> None:
        """Run the annealing / MCMC loop."""
        # Book‑keeping
        self.energies_: Deque[float] = deque()
        self.acceptance_count_ = 0
        self.total_proposals_ = 0
        self.current_graph_ = self.chain.state
        self.current_energy_ = self.energy_fn(self.current_graph_)

        # Initialize best graph along with it's energy
        self.best_graph_ = (self.current_graph_, self.current_energy_)

        for step in tqdm(range(self.max_steps)):
            beta = self.schedule(step)
            new_graph = self.chain._propose()
            new_energy = self.energy_fn(new_graph)

            # Metropolis acceptance probability
            dE = new_energy - self.current_energy_
            accept = dE < 0 or self.random.random() < exp(-beta * dE)

            self.total_proposals_ += 1
            if accept:
                self.acceptance_count_ += 1
                self.current_graph_, self.current_energy_ = new_graph, new_energy
                self.chain.state = self.current_graph_

                # Update best graph if new energy minimum is achieved
                if self.best_graph_[1] > self.current_energy_:
                    self.best_graph_ = (self.current_graph_, self.current_energy_)

            self.energies_.append(self.current_energy_)