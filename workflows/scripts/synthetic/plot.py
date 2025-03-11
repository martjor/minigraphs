import matplotlib.pyplot as plt 
from scripts.utils.io import load_dict
import numpy as np

data = {key:[] for key in ['er','ws','ba']}
for file in snakemake.input:
    model, n_nodes = file.split('/')[-2].split('_')
    n_nodes = int(n_nodes)
    eig_1 = load_dict(file)['eig_1']
    
    # Store in dictionary
    data[model].append((n_nodes, eig_1))
    
fig, ax = plt.subplots(figsize=(10,5),dpi=300)
for model, data in data.items():
    data = np.array(data)
    
    ax.scatter(
        data[:,0],
        data[:,1],
        label=model
    )
    
ax.legend()
ax.set_title("Spectral Radii of Synthetic Graphs")
ax.set_xlabel("Number of Nodes")
ax.set_ylabel("Dominant Eigenvalue")
    
plt.savefig(snakemake.output[0])
    