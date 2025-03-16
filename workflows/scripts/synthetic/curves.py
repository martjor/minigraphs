import matplotlib.pyplot as plt 
import numpy as np

# Collect inputs into dictionary
data = {key:{'data':[],'n_nodes':[]} for key in ['er','ws','ba']}

for file in snakemake.input:
    model, n_nodes = file.split('/')[-2].split('_')
    n_nodes = int(n_nodes)
    
    # Load infected data for graph
    data[model]['data'].append(
        np.load(file)[:,1,:].mean(0)
    )
    
    # Save number of nodes
    data[model]['n_nodes']
    
fig, axes = plt.subplots(1,4,figsize=(20,5),dpi=300)

for i, model in enumerate(data.keys()):
    for j, trajectory in enumerate(data[model]['data']):
        axes[i].plot(
            trajectory,
            c='red',
            alpha=0.5
        )    
        
    axes[i].set_title(f"{model.capitalize()}")
    axes[i].set_ylim([0,1])
    axes[i].set_xlim([0,150])
        
plt.savefig(snakemake.output[0])

    