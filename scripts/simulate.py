import networkx as nx 
from numpy.random import rand

class Contagion:
    def __init__(self,graph,
                      node_r,
                      node_b,
                      n_iterations):
        self.graph = graph 
        self.node_r = node_r 
        self.node_b = node_b 
        self.n_iterations = n_iterations
        
        # Generalized adoption function
        self.func_h = lambda a_r, a_b: a_r and (a_r / (a_r + a_b))
    
    def func_selection(self,a_r,a_b):
        '''Selection function
        '''
        probability = self.func_h(a_r,a_b) / (self.func_h(a_r,a_b) +
                                              self.func_h(a_b,a_r))
        
        return probability 
    
    def func_switching(self,a_r,a_b):
        '''Switching function
        '''
        probability = self.func_h(a_r,a_b) + self.func_h(a_b,a_r)
        
        return probability
    
    def infect(self,node,color):
        '''Infects a node to the specified color and increases the count
        '''
        nx.set_node_attributes(self.graph,{node: {'color':color}})
        self.count[color] += 1
        
    def update(self,node):
        '''Updates the state of a node based on its neighborhood
        '''
        if self.graph.nodes[node]['color'] == 'gray':
            alpha_r = 0
            alpha_b = 0
            
            neighborhood = self.graph[node]
            for neighbor in neighborhood:
                if self.graph.nodes[neighbor]['color'] == 'red':
                    alpha_r += 1
                elif self.graph.nodes[neighbor]['color'] == 'blue':
                    alpha_b += 1 
                    
            if rand() < self.func_switching(alpha_r,alpha_b):
                if rand() < self.func_selection(alpha_r,alpha_b):
                    self.infect(node,'red')
                else:
                    self.infect(node,'blue')
            
    
    def run(self):
        '''Simulates contagion on the graph
        '''
        self.count = {'red':0,'blue':0}
        
        # Initialize graph colors
        nx.set_node_attributes(self.graph,'gray','color')
            
        # Initial infection
        for node in set([self.node_r, self.node_b]):            
            if self.node_r == self.node_b:
                if rand() < 0.5:
                    self.infect(node,'blue')
                else: 
                    self.infect(node,'red')
            elif self.node_r == node:
                self.infect(node,'red')
            else:
                self.infect(node,'blue')
                    
        
        # Simulate infection
        iteration = 0

        while iteration < self.n_iterations:
            for node in self.graph.nodes:
                # Update node's state
                self.update(node)
                
            iteration += 1
                

            
        
        
        
            
    
            
        
    
    
    




