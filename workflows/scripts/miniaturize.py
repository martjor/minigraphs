from mpi4py.MPI import COMM_WORLD
from minigraphs.mcmc import ParallelTempering
from minigraphs.mcmc.chains import SubgraphBoundary
from minigraphs.metrics import spectral_radius
import networkx as nx 
from minigraphs.data import load_graph
import click

@click.command
@click.argument('betas', type=click.FLOAT, nargs=-1)
@click.option('--n-steps', type=click.INT)
@click.option('--exchange-freq', type=click.INT)
@click.option('--n-nodes', type=click.INT)
@click.option('--out-graph', type=click.Path())
@click.option('--out-df', type=click.Path())
def main(
        betas,
        n_steps, 
        exchange_freq,
        n_nodes,
        out_graph, 
        out_df,
    ):
    n_swaps     = 10
    hamsterster = load_graph('hamsterster')

    # Initialize annealer
    runner = ParallelTempering(
        COMM_WORLD,
        annealer_data=[(SubgraphBoundary(hamsterster, n_nodes, n_swaps), inv_temp) for inv_temp in betas],
        energy=lambda graph: -spectral_radius(graph),
        exchange_freq=exchange_freq,
        n_steps=n_steps,
    )

    # Miniaturize
    runner.run()

    # Save graph
    runner.best_graph_save(nx.write_adjlist, {'path': out_graph})

    # Save miniaturization data
    runner.gather_results()
    runner.results_to_csv(out_df)
    
if __name__ == '__main__':
    main()
