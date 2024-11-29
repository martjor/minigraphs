import numpy as np
import matplotlib.pyplot as plt

matrix = np.load(snakemake.input[0])

fig, ax = plt.subplots(figsize=(5,5),dpi=300)
im = ax.imshow(matrix)
fig.tight_layout()

plt.savefig(snakemake.output[0])