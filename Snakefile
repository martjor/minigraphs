configfile: "config/config.yaml"
import numpy as np
import yaml
import glob
import pandas as pd


graph_name = config['name']
n_vertices = config['generator']['n_vertices']
rule all:
    input:
        "results/test/equilibria.csv"

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

checkpoint find_equilibria:
    input:
        "data/{name}/payoff_red.npy",
        "data/{name}/payoff_blue.npy"
    output:
        directory("data/{name}/equilibria")
    script:
        "scripts/equilibria.py"

rule simulate_equilibrium:
    input:
        "{path}/adjacency.npz",
        "{path}/equilibria/equilibrium_{idx}.npy"
    output:
        "{path}/equilibria/count_{idx}.yaml",
        "{path}/equilibria/graph_{idx}.gexf"
    params:
        n_iterations=config['payoff_equilibria']['n_iterations'],
        n_trials=config['payoff_equilibria']['n_trials']
    script:
        "scripts/compute_payoff_equilibria.py"

def equilibria(wildcards):
    equilibria_dir = checkpoints.find_equilibria.get(**wildcards).output[0]

    idx, = glob_wildcards(f"data/{wildcards.name}/equilibria/equilibrium_{{idx}}.npy")
    files = expand(f"data/{wildcards.name}/equilibria/count_{{idx}}.yaml",idx=idx)
    return files

rule aggregate_equilibria:
    input:
        equilibria 
    output:
        "results/{name}/equilibria.csv"
    params:
        indices=lambda w: glob_wildcards(f"data/{w.name}/equilibria/equilibrium_{{idx}}.npy").idx
    run:
        dicts = []
        for file in input:
            with open(file,'r') as f:
                dicts.append(yaml.safe_load(f))
        
        df = pd.DataFrame(dicts)
        df['idx'] = params.indices
        df.set_index('idx',inplace=True)

        df['welfare'] = df['red'] + df['blue']

        df.to_csv(output[0])

rule price_of_anarchy:
    input:
        red="data/{name}/payoff_red.npy",
        blue="data/{name}/payoff_blue.npy",
        equilibria="results/{name}/equilibria.csv"
    notebook:
        "notebooks/analysis.py.ipynb"






