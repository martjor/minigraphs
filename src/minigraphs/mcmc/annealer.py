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

        # Book‑keeping
        self.samples: Deque[Graph] = deque()
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
            accept = dE < 0 or self.random.random() < exp(-beta * dE)

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
    def __iter__(self) -> Iterable[Graph]:
        return iter(self.samples)
