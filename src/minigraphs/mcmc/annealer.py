from typing import Union, Callable, Optional
from collections import deque
from tqdm import tqdm 
from .chains import Chain 
from math import exp
import random 
from networkx import Graph
from pandas import DataFrame
from typing import NamedTuple

class State(NamedTuple):
    """A named tuple holding the current state of the annealer.
    """
    beta: float 
    energy: float

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
    verbose : bool, default=False

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
        n_steps: int = 1_000,
        seed: Optional[int] = None,
        verbose: Optional[bool]= False,

    ) -> None:
        self.chain = chain
        self.energy_fn = energy
        self.schedule = schedule 
        self.n_steps = n_steps
        self.random = random.Random(seed)
        self.verbose = verbose

    def run(self) -> None:
        """Run the annealing / MCMC loop."""
        self._history_ = deque()

        # Initialize internal state
        old_graph = self.chain.state 
        old_energy = self.energy_fn(old_graph)

        # Initialize schedule
        schedule = self.schedule if callable(self.schedule) else lambda _: float(self.schedule)

        # Initialize best graph along with it's energy
        self.best_graph_ = old_graph
        self.best_energy_ = old_energy

        for step in tqdm(range(self.n_steps), disable=not self.verbose):
            beta = schedule(step)
            new_graph = self.chain._propose()
            new_energy = self.energy_fn(new_graph)

            # Metropolis acceptance probability
            dE = new_energy - old_energy
            accept = dE < 0 or self.random.random() < exp(-beta * dE)

            if accept:
                # Update chain
                self.chain.state = new_graph

                # Update internal variables
                old_graph, old_energy = new_graph, new_energy

                # Update best graph if new energy minimum is achieved
                if self.best_energy_ > new_energy:
                    self.best_graph_ = new_graph
                    self.best_energy_ = new_energy

            # Store state
            self._history_.append(
                State(
                    beta=beta,
                    energy=old_energy
                )
            )

    @property
    def history_(self) -> DataFrame:
        """Retrieves the annealer's state history.
        """
        return DataFrame(self._history_,)