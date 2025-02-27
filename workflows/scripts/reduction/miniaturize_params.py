#!/usr/bin/env python3
from minigraphs.miniaturize import MH
from minigraphs.graph import spectral_radius

from numpy import log, inf, nan
import pandas as pd 
import networkx as nx
import sys
import os
import yaml
from scripts.utils.io import StreamToLogger
import logging
from scripts.reduction.pt_setup import DICT_METRICS_FUNCS
'''Calculates the parameters for the specified graph
'''

# Get the log file from Snakemake
log_file = snakemake.log[0]

# Configure logging to write to the Snakemake log file
logging.basicConfig(
    filename=log_file,
    level=logging.DEBUG,  # Capture all logs (DEBUG and above)
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Replace stdout and stderr with the logger
sys.stdout = StreamToLogger(logging.getLogger(), logging.INFO)
sys.stderr = StreamToLogger(logging.getLogger(), logging.ERROR)

def weights(
    metrics_file,
    params_file,
    n_vertices,
    n_iterations,
    n_trials,
    targets
    ):
    
    # Retrieve Graph Metrics
    with open(metrics_file) as file:
        graph_metrics = yaml.safe_load(file)
        
    # Construct miniaturization target metrics and functions
    metrics={}
    functions={}
    for target in targets:
        metrics[target] = graph_metrics[target]
        functions[target] = DICT_METRICS_FUNCS[target]
    
    print(f"Calculating parameters for graph at graph at {metrics_file}")
    print(f"\t - Size: {n_vertices} nodes ({1 - n_vertices/graph_metrics["n_nodes"]:.02f}% miniaturization)")
    print(f"\t - Number iterations per sample: {n_iterations}")
    print(f"\t - Number of samples: {n_trials}\n")

    # Calculate weights
    params = []
    for i in range(n_trials):
        print(f"Sweep {i+1}/{n_trials}")
        
        # Construct replica
        replica = MH(functions,
                     schedule=lambda beta:0,
                     n_changes=1)
        
        # Transform ER graph
        G = nx.erdos_renyi_graph(n_vertices,graph_metrics['density'])
        replica.transform(G,
                          metrics,
                          n_iterations=n_iterations)

        # Retrieve trajectories
        df = replica.trajectories_
        weights = dict(1/df[replica.metrics].diff().abs().mean())
        
        print("Weights:")
        print(weights)

        # Calculate optimal beta
        replica = MH(functions,
                     schedule=lambda beta:0,
                     weights=weights)
        
        # Transform ER graph
        G = nx.erdos_renyi_graph(n_vertices,graph_metrics['density'])
        replica.transform(G,
                          metrics,
                          n_iterations=n_iterations)

        df = replica.trajectories_

        beta = -log(0.23) * 1/df['Energy'].diff().abs().mean()

        print(f"Beta: {beta}\n")

        params.append([beta] + list(weights.values()))

    params = pd.DataFrame(params,columns=['beta'] + list(weights.keys()))
    print("Measeured parameters:")
    print(params,f"\n")

    params.replace([inf, -inf], nan, inplace=True)
    params = params.mean()
    print("Final parameters:")
    print(params)

    # Save to yaml
    params = params.to_dict()
    params_dict = {'beta': params['beta']}
    params.pop('beta')
    params_dict['weights'] = params 
    
    with open(params_file,'w') as file:
        yaml.dump(params_dict,file,default_flow_style=False)
        
weights(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.params.n_nodes,
    snakemake.params.n_steps, 
    snakemake.params.n_trials,
    snakemake.params.targets,
)
    