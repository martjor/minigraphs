from minigraphs import simulation as sim
from scipy.sparse import load_npz
import numpy as np
import yaml
 
# Parameters
params = snakemake.params.sir_params

sir = sim.Sir(params['tau'], params['gamma'])

# Instantiate simulation object
simulation = sim.Simulation(load_npz(snakemake.input[0]))

# Allocate memory
shape = (params['n_trials'], 3, params['n_steps'])
results = np.zeros(shape)

# Simulate epidemic
for i in range(params['n_trials']):
    simulation.run(sir, params['n_steps'])
    
    results[i,:,:] = simulation.trajectories_.T
    
# Save simulation results
np.save(snakemake.output[0],results)