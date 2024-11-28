configfile: "config/config.yaml"
import numpy as np

graph_name = config['name']

rule all:
    input:
        "data/test/adjacency.npz"

rule generate_graph:
    output:
        "data/{name}/adjacency.npz"
    params:
        params_gen=config['generator']
    script:
        "scripts/generate.py"

# rule payoff:
#     input:
#         "{path}/adjacency.npz"
#     output:
#         temp("{path}/entry_{i}_{j}.npy")
#     params:
#         n_trials=config['payoff']['n_trials']
#     script:
#         "scripts/simulate.py"

# rule payoff:
#     input:
#         entries=[f"data/entry_{i}_{j}.npy" for i in range(n_vertices) for j in range(n_vertices)]
#     output:
#         "data/payoff.npy"
#     run:
#         # Allocate memory for payoff matrix
#         M = np.zeros((n_vertices,n_vertices))

#         # Read each entry into the matrix
#         for i in range(n_vertices):
#             for j in range(n_vertices):
#                 idx = i * n_vertices + j
#                 M[i,j] = np.load(input.entries[idx])

#         # Save matrix
#         np.save(output[0],M)

# rule draw:
#     input:
#         "data/payoff.npy"
#     output:
#         "results/payoff.png"
#     script: 
#         "scripts/draw_payoff.py"


