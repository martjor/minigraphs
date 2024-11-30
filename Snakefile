configfile: "config/config.yaml"
import numpy as np
import yaml


graph_name = config['name']
n_vertices = config['generator']['n_vertices']
rule all:
    input:
        "results/test/payoffs.png",
        "data/test/equilibria/"

rule generate_graph:
    output:
        "data/{name}/adjacency.npz"
    params:
        params_gen=config['generator']
    script:
        "scripts/generate.py"

rule payoff_entries:
    input:
        "{path}/adjacency.npz"
    output:
        "{path}/entries/count_{i}_{j}.yaml",
        "{path}/entries/graph_{i}_{j}.gexf"
    
    params:
        n_trials=config['payoff']['n_trials'],
        n_iterations=config['payoff']['n_iterations']
    script:
        "scripts/compute_payoff.py"

rule payoff_matrix:
    input:
        entries=lambda w:[f"data/{w.name}/entries/count_{i}_{j}.yaml" for i in range(n_vertices) for j in range(n_vertices)]
    output:
        "data/{name}/payoff_red.npy",
        "data/{name}/payoff_blue.npy"
    run:
        # Allocate memory for payoff matrix
        red = np.zeros((n_vertices,n_vertices))
        blue = np.zeros((n_vertices,n_vertices))

        # Read each entry into the matrix
        for i in range(n_vertices):
            for j in range(n_vertices):
                idx = i * n_vertices + j

                with open(input.entries[idx],'r') as file:
                    count = yaml.safe_load(file)

                red[i,j] = count['red']
                blue[i,j] = count['blue']

        # Save matrix
        np.save(output[0],red)
        np.save(output[1],blue)

rule payoff_draw:
    input:
        "data/{name}/payoff_red.npy",
        "data/{name}/payoff_blue.npy"
    output:
        "results/{name}/payoffs.png"
    script: 
        "scripts/payoff_draw.py"

rule find_equilibria:
    input:
        "{path}/payoff_red.npy",
        "{path}/payoff_blue.npy"
    output:
        directory("{path}/equilibria/")
    script:
        "scripts/equilibria.py"

# rule simulate_equilibrium:
#     input:
#         "{path}/equilibrium_{idx}.npy"
#     output:
#         "{path}/count_{idx}.yaml",
#         "{path}/graph_{idx}.gexf"
#     script:
#         "scripts/"






