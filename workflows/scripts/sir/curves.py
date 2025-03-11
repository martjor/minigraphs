import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np


fig, ax = plt.subplots(figsize=(8,4), dpi=300)
cmap = plt.get_cmap('Reds')

# Plot Original Curve
file = snakemake.input.original

infected = np.load(snakemake.input.original)[:,1,:].mean(0)
peak_original = np.array([infected.argmax(), infected.max()])

arrowprops={
    'arrowstyle':"->",
    'connectionstyle':'arc3'
}

ax.plot(
    infected,
    color=cmap(1.0),
    label="original",
)

xytext = peak_original * [10.0, 1.0]
ax.annotate(
    r'$I_{\max}$',
    xy=peak_original,
    xytext=xytext,
    arrowprops=arrowprops
)


data = {}
for file in snakemake.input.miniatures:
    miniature = '_'.join(file.split('/')[-2].split('_')[:-1])
    infected = np.load(file)[:,1,:]
    
    if data.get(miniature) is None:
        data[miniature] = [infected,]
    else:
        data[miniature].append(infected)

sizes = dict(zip(data.keys(),[int(key.split('_')[1]) for key in data.keys()]))
sizes_list = sorted(sizes.values(), reverse=True)
max_nodes = max(sizes.values())
min_nodes = min(sizes.values())

for file in sorted(snakemake.input.ers):
    n_nodes = int(file.split('/')[-2].split('_')[1])
    color = 0.2 + 0.6 * (n_nodes-min_nodes) / (max_nodes-min_nodes)
    infected = np.load(file)[:,1,:].mean(0)
    
    x = np.arange(len(infected)) + 100 * (sizes_list.index(n_nodes) + 1)
    ax.plot(
        x,
        infected,
        linestyle='--',
        color=cmap(color),
        alpha=0.5
    )

for miniature, list in sorted(data.items(),reverse=True):
    color = 0.2 + 0.6 * (sizes[miniature]-min_nodes) / (max_nodes-min_nodes)
    
    infected = np.array(list).mean((0,1))
    peak = np.array([infected.argmax(), infected.max()])

    if peak[1] > 0.05:
        ax.annotate(
        f"{peak[1]/peak_original[1]:.2f} " + r"$I_{\max}$", 
        xy=peak,
        xytext=(xytext[0], peak[1]),
        arrowprops=arrowprops
    )
    
    ax.plot(
        infected,
        label=sizes[miniature],
        color=cmap(color)
    )

    
ax.legend()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylabel("Normalized Population"),
ax.set_xlabel("Time")
ax.set_title(f"Infected curve on {snakemake.wildcards.network} and miniatures")



plt.savefig(snakemake.output[0])