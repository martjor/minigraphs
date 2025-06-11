import mesa
from mesa.space import NetworkGrid
from networkx import Graph
from enum import Enum

class Compartment(Enum):
    S = 0
    I = 1
    R = 2

class IndividualAgent(mesa.Agent):
    def __init__(self, model: "SIRModel", compartment: Compartment, node: int):
        super().__init__(model)
        self.compartment = compartment
        self._compartment = None 
        self.node = node
        self.got_infected = False

    def compute_state(self):
        '''Computes what the next state of the agent should be
        '''
        match self.compartment:
            case Compartment.S:
                # Get neighbors of agent
                neighbors = self.model.network.get_neighbors(self.node)

                # Count infected neighbors
                n_infected = 0 
                for neighbor in neighbors:
                    n_infected += 1 if neighbor.compartment == Compartment.I else 0

                # Update state according to probability 
                p_infected = 1 - (1 - self.model.beta) ** n_infected

                if self.model.random.random() < p_infected:
                    self._compartment = Compartment.I

                    # Update flag
                    self.got_infected = True
                else:
                    self._compartment = Compartment.S

                
            
            case Compartment.I:
                # Update flag
                if self.got_infected:
                    self.got_infected = False 

                # Calculate probability of recovery
                if self.model.random.random() < self.model.gamma:
                    self._compartment = Compartment.R

            case Compartment.R:
                # Recovered individuals remain recovered
                pass

    def update_state(self):
        '''Updates state to the previously computed state
        '''
        self.compartment = self._compartment 

def count_compartment(model: "SIRModel", compartment: Compartment) -> int:
    count = 0
    for agent in model.agents:
        count += 1 if agent.compartment == compartment else 0

    return count

def count_new_infections(model: "SIRModel") -> int:
    count = 0
    for agent in model.agents:
        count += agent.got_infected

    return count
        

class SIRModel(mesa.Model):
    """Agent-based SIR model for the spread of disease.

    Parameters
    ----------
    beta: float
        Infection probability over a link.
    gamma: float
        Recovery rate of an infected individual.
    network: nx.Graph
        Underlying contact network.
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

        # Create agents
        compartments = [Compartment.S] * (self.n_agents - n_infected) + [Compartment.I] * n_infected
        self.random.shuffle(compartments)

        agents = IndividualAgent.create_agents(
            self,
            n=self.n_agents,
            compartment=compartments,
            node=list(network.nodes())
        )

        # Create network and add agents to network
        self.network = NetworkGrid(network)
        for agent in agents: 
            self.network.place_agent(agent, agent.node)

        self.datacollector = mesa.DataCollector(
            model_reporters={
                'S': lambda model: count_compartment(model, Compartment.S),
                'I': lambda model: count_compartment(model, Compartment.I),
                'R': lambda model: count_compartment(model, Compartment.R),
                'n_new_infections': lambda model: count_new_infections(model)
            }
        )

    def step(self):
        """Advances the model one time step.
        """ 
        self.datacollector.collect(self)
        self.agents.do('compute_state') 
        self.agents.do('update_state')