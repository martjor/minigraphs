from mpi4py import MPI
from .annealer import SimulatedAnnealing
from .chains import Chain
from typing import List, Tuple, Callable
from networkx import Graph
from math import ceil, exp
from collections import deque
import pandas as pd
import random


class ParallelTempering:
    def __init__(
            self,
            annealer_data: List[Tuple[Chain, float]],
            energy: Callable[[Graph], float],
            exchange_freq: int, 
            n_steps: int, 
    ):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        self.annealer_data = annealer_data
        self.energy = energy
        self.exchange_freq = exchange_freq
        self.n_steps = n_steps

    def _initialize(self):
        """Initializes the annealer for the current replica.
        """

        # Validate annealer data
        if (len(self.annealer_data) != self.size) and (self.rank == 0):
            raise ValueError(f"Number of annealers provided ({len(self.annealer_data)}) does not match size of comm ({self.size})")
        
        chain, self.inv_temp = sorted(self.annealer_data, key=lambda tuple: tuple[1])[self.rank]

        # Initialize annealer
        self.annealer = SimulatedAnnealing(
            chain,
            self.energy,
            self.inv_temp
        )

        self.history_ = deque()

    def _attempt_exchange(self, my_energy: float, my_inv_temp: float) -> float:
        """Replicas exchange inverse temperatures probabilistically.
        """
        partner = self.rank + 1 if self.rank % 2 == 0 else self.rank - 1

        if 0 <= partner < self.size:
            # Exchange energy and inverse temperatures. 
            partner_energy, partner_inv_temp = self.comm.sendrecv(
                (my_energy, my_inv_temp),
                dest=partner, 
                source=partner,
            )

            dB = partner_inv_temp - my_inv_temp
            dE = partner_energy - my_energy

            if dE < 0 or random.random() < exp(-dB*dE):
                return partner_inv_temp
            else:
                return my_inv_temp
            
    def run(self):
        """Runs the parallel tempering algorithm.
        """
        # Initialize properties of the annealer.
        self._initialize()

        n_episodes = ceil(self.n_steps / self.exchange_freq)
        steps = 0

        for episode in range(n_episodes):
            # Update number of steps per episode
            n_steps_per_episode = self.exchange_freq if episode < (n_episodes-1) else (self.n_steps % self.exchange_freq)

            # Run annealer
            self.annealer.max_steps = n_steps_per_episode
            self.annealer.run()
            steps += n_steps_per_episode

            # Append history
            self.history_.append(
                (n_steps_per_episode, self.inv_temp, list(self.annealer.energies_))
            )

            # Attempt exchange
            self.inv_temp = self._attempt_exchange(self.annealer.current_energy_, self.inv_temp)
            self.annealer.schedule = self.inv_temp

    def gather_results(self):
        """Gathers the histories from all replicas.
        """
        history = list(self.history_)
        all_histories = self.comm.gather(history, root=0)

        if self.rank == 0:
            df_list = []

            for rid, hist in enumerate(all_histories):
                df = pd.DataFrame(hist, columns=["step", "temp", "energy"])
                df["replica"] = rid 
                df_list.append(df)

            return pd.concat(df_list, ignore_index=True)
        return None

    




