from mpi4py.MPI import COMM_WORLD
from minigraphs.mcmc import ParallelTempering
from minigraphs.mcmc.chains import SubgraphBoundary
from minigraphs.metrics import graph_spectrum
import networkx as nx 
from minigraphs.data import load_graph
from seaborn import lineplot
from matplotlib.pyplot import figure
import numpy as np 

hamsterster = load_graph('hamsterster')
inverse_temperatures = [10, 100, 300, 500, 700, 1000]
n_nodes = 100

def spectral_radius(graph: nx.Graph) -> float:
    """Calculates the spectral radius of a graph."""
    return graph_spectrum(graph)[0]

runner = ParallelTempering(
    COMM_WORLD,
    annealer_data=[(SubgraphBoundary(hamsterster, n_nodes, 1), inv_temp) for inv_temp in inverse_temperatures],
    energy=lambda graph: -spectral_radius(graph),
    exchange_freq=100,
    n_steps=30000,
)

runner.run()
runner.gather_results()
runner.results_to_csv("test.csv")

if runner.rank == 0:
    history = runner._results.reset_index()

    figure()
    g = lineplot(
        history,
        x='step',
        hue='replica',
        y='energy'
    )

    g.figure.savefig('energy.png')

    figure()
    g = lineplot(
        history,
        x='step',
        hue='replica',
        y='beta'
    )
    
    g.figure.savefig('schedule.png')

runner.best_graph_save(nx.write_adjlist, path='best_graph.adjlist')