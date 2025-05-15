import mesa
from mesa.space import NetworkGrid
from networkx import Graph

class IndividualAgent(mesa.Agent):
    def __init__(self, model, compartment, node):
        super().__init__(model)
        self.compartment = compartment
        self._compartment = None 
        self.node = node

    def compute_state(self):
        '''Computes what the next state of the agent should be
        '''
        match self.compartment:
            case "S":
                # Get neighbors of agent
                neighbors = self.model.network.get_neighbors(self.node)

                # Count infected neighbors
                n_infected = 0 
                for neighbor in neighbors:
                    n_infected += 1 if neighbor.compartment == "I" else 0

                # Update state according to probability 
                p_infected = 1 - (1 - self.model.beta) ** n_infected

                if self.model.random.random() < p_infected:
                    self._compartment = "I"

                    self.model.compartment_count["S"] -= 1
                    self.model.compartment_count["I"] += 1
                else:
                    self._compartment = "S" 
            
            case "I":
                # Calculate probability of recovery
                if self.model.random.random() < self.model.gamma:
                    self._compartment = "R" 

                    self.model.compartment_count["I"] -= 1
                    self.model.compartment_count["R"] += 1

            case _:
                # Recovered individuals remain recovered
                pass

    def update_state(self):
        '''Updates state to the previously computed state
        '''
        self.compartment = self._compartment 


class SIRModel(mesa.Model):
    """Agent-based SIR model for the spread of disease.

    Parameters
    ----------
    beta: float
        Infection probability over a link.
    gamma: float
        Recovery rate of an infected individual.
    network: nx.Graph
        Underlying interaction network.
    n_infected: int
        Number of infected individuals at the beginning of the simulation.
    seed: float
        Random state of the model
    """
    def __init__(
            self, 
            beta: float, 
            gamma: float, 
            network: Graph, 
            n_infected: int=1, 
            seed: int=None
        ):
        super().__init__(seed=seed)

        self.n_agents = network.number_of_nodes()
        self.beta = beta 
        self.gamma = gamma

        # Initialize compartment count
        self.compartment_count = dict.fromkeys(["S","I","R"], 0)
        self.compartment_count["S"] = self.n_agents - n_infected
        self.compartment_count["I"] = n_infected

        # Create agents
        compartments = ["S"] * (self.n_agents - n_infected) + ["I"] * n_infected
        self.random.shuffle(compartments)

        agents = IndividualAgent.create_agents(
            self,
            n=self.n_agents,
            compartment=compartments,
            node=list(range(self.n_agents))
        )

        # Create network and add agents to network
        self.network = NetworkGrid(network)
        for agent in agents: 
            self.network.place_agent(agent, agent.node)

        self.datacollector = mesa.DataCollector(
            model_reporters={
                'S': lambda m: m.compartment_count['S'],
                'I': lambda m: m.compartment_count['I'],
                'R': lambda m: m.compartment_count['R'],
            }
        )

    def step(self):
        """Advances the model one time step.
        """ 
        self.datacollector.collect(self)
        self.agents.do('compute_state') 
        self.agents.do('update_state')