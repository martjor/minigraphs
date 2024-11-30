import nashpy as nash 
import numpy as np 
import os

red = np.load(snakemake.input[0])
blue = np.load(snakemake.input[1])

game = nash.Game(red,blue)
os.makedirs(snakemake.output[0],exist_ok=True)

for i, eq in enumerate(game.vertex_enumeration()):
    # Get rid of negative numbers
    eq = np.abs(eq)

    # Store nash equilibrium
    np.save(os.path.join(snakemake.output[0],f"equilibrium_{i}.npy"),eq)

