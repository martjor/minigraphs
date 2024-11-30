import matplotlib 
import matplotlib.pyplot as plt 
import numpy as np 

red = np.load(snakemake.input[0])
blue = np.load(snakemake.input[1])

fig, axes = plt.subplots(1,2,dpi=300,figsize=(10,5))
axes[0].imshow(red)
axes[0].set_title("Red")

axes[1].imshow(blue)
axes[1].set_title("Blue")

plt.savefig(snakemake.output[0],bbox_inches='tight')
