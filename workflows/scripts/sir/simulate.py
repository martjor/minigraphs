from minigraphs import simulation as sim
from scipy.sparse import load_npz
import numpy as np
import yaml
 
# Parameters
tau = snakemake.params.tau
gamma = snakemake.params.gamma
n_steps = snakemake.params.n_steps
n_trials = snakemake.params.n_trials

sir = sim.Sir(tau, gamma)

# Instantiate simulation object
simulation = sim.Simulation(load_npz(snakemake.input[0]))

# Allocate memory
shape = (n_trials, 3, n_steps)
results = np.zeros(shape)

# Simulate epidemic
for i in range(n_trials):
    simulation.run(sir, n_steps)
    
    results[i,:,:] = simulation.trajectories_.T
    
# Save simulation results
np.save(snakemake.output[0],results)