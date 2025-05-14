from mpi4py.MPI import Comm
from .annealer import SimulatedAnnealing
from .chains import Chain
from typing import List, Tuple, Callable, Optional
from networkx import Graph
from math import ceil, exp
from collections import deque
import pandas as pd
import random
from functools import wraps
from tqdm import tqdm 

class ParallelTempering:
    def __init__(
            self,
            comm: Comm,
            annealer_data: List[Tuple[Chain, float]],
            energy: Callable[[Graph], float],
            exchange_freq: int, 
            n_steps: int,
            *, 
            verbose: Optional[bool]=False,
    ):
        """Implementation of the Parallel Tempering algorithm on graphs.

        Various replicas (instances of the `SimulatedAnnealing`) class are instantiated and allowed to optimize
        the graph structure. These replicas will exchange inverse temperatures every `exchange_freq` steps.

        Parameters
        ----------
        comm : Comm 
            MPI group coordinating each rank.
        annealer_data : List[Tuple[Chain, float]]
            List of tuples containing the chain to be used by each replica, along with their corresponding inverse temperatures.
        energy : Callable[[Graph], float]
            Function that calculates the energy of a graph.
        exchange_freq : int
            Number of steps between exchanges.
        n_steps:
            Total number of steps to run the algorithm.
        """
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        self.annealer_data = annealer_data
        self.energy = energy
        self.exchange_freq = exchange_freq
        self.n_steps = n_steps
        self.verbose = verbose

    def _exchange(self, partner: int, my_inv_temp: float) -> float:
        """Proposes swap to partner by sending them current energy and inv temp.

        Replicas send their inverse temperature if the rank of their partner is higher (`partner > self.rank`)
        and receive otherwise:
            - Receiver (higher rank) decides whether to swap.
            - Sender (lower rank) adopts the returned inverse temperature unconditionally.

        Parameters
        ----------
        partner : int
            Partner to exchange inverse temperatures with.
        my_inv_temp : float 
            Current inverse temperature ofo the replica.

        Returns
        -------
        float
            The new inverse temperature to use.
        """
        my_energy = self.annealer._history_[-1].energy

        # Send information
        if (partner == (self.rank+1)) and (partner < self.size):
            # Send energy and inverse temperature for evaluation
            self.comm.send(
                (my_energy, my_inv_temp),
                dest=partner, 
            )

            # Receive new inverse temperature from partner
            my_inv_temp = self.comm.recv(
                source=partner,
            )
            

        # Receive information from partner
        elif (partner == (self.rank-1) and (partner >= 0)):
            partner_energy, partner_inv_temp = self.comm.recv(
                source=partner,
            )

            dB = partner_inv_temp - my_inv_temp
            dE = partner_energy - my_energy
            exponent = dB * dE

            # Accept or reject inv temperature swap
            accept = exponent < 0 or random.random() < exp(-exponent)
            if accept:
                my_inv_temp, partner_inv_temp = partner_inv_temp, my_inv_temp

            # Send inverse temperature back to partner
            self.comm.send(
                partner_inv_temp,
                dest=partner,
            )

        return my_inv_temp
            
    def run(self):
        """Runs the parallel tempering algorithm.
        """
        # Initialize properties of the annealer.
        if len(self.annealer_data) != self.size:
            raise ValueError(f"Mismatch between annealer data ({len(self.annealer_data)}) and comm size ({self.size})")

        sorted_data = sorted(self.annealer_data, key=lambda item: item[1])
        chain, beta = sorted_data[self.rank]

        # Initialize annealer
        self.annealer = SimulatedAnnealing(
            chain,
            self.energy,
            beta
        )

        # Estimate number of episodes
        n_episodes = ceil(self.n_steps / self.exchange_freq)
        n_steps_remain = self.n_steps

        # Allocate memory for history
        self._history_ : List[pd.DataFrame] = [None] * n_episodes

        for episode in tqdm(range(n_episodes)):
            # Number of steps in the current episode
            remainder = n_steps_remain % self.exchange_freq
            n_steps = self.exchange_freq if remainder == 0 else remainder
            
            # Run annealer
            self.annealer.n_steps = n_steps
            self.annealer.run()
            n_steps_remain -= n_steps

            # Append history
            self._history_[episode] = self.annealer.history_
            
            # Even ranks initiate exchange
            partner = self.rank + 1 if (self.rank % 2 == 0) else self.rank - 1
            beta = self._exchange(partner, beta)

            # Odd ranks initiate exchange
            partner = self.rank + 1 if (self.rank % 2 == 1) else self.rank - 1
            beta = self._exchange(partner, beta)

            self.annealer.schedule = beta
                
        self.comm.Barrier()        

    @property
    def history_(self) -> pd.DataFrame:
        """History of the replica.
        """
        return pd.concat(self._history_)

    def gather_results(self):
        """
        Gathers the histories from all replicas and stores them into a `pd.Series` object.
        """
        if (self.rank == 0) and self.verbose:
            print("Collecting histories...")

        histories = self.comm.gather(
            self.history_.to_dict('list'), 
            root=0
        )

        if self.rank == 0:
            self._results = pd.concat(
                (pd.DataFrame(history) for history in histories),
                keys=range(len(histories)),
                names=['replica','step']
            )

    @wraps(pd.DataFrame.to_csv)
    def results_to_csv(self, *args, **kwargs):
        """Saves the trajectories of each replica to a csv file.
        """
        if self.rank == 0:
            self._results.to_csv(*args, **kwargs)

    

    




