import pandas as pd
import networkx as nx
from typing import Dict, Callable
from collections import deque
from tqdm import tqdm 
from .mcmc.chains import Chain

def subgraph_metrics(
        chain: Chain,
        n_samples: int,
        metrics: Dict[str, Callable[[nx.Graph], float]]
    ) -> pd.DataFrame:
    # Empty list for metrics
    measurements = deque()

    for _ in tqdm(range(n_samples)):
        graph = next(chain)

        measurements.append(
            {metric: func(graph) for metric, func in metrics.items()}
        )

    return pd.DataFrame(measurements)