from multiprocessing import Pool
import networkx as nx
from itertools import combinations
import string
import numpy as np
import random
import json
import pickle
import copy


#model constructor
class SimplagionModel():
    def __init__(self, node_neighbors_dict, triangles_list, I_percentage):
        
        #parameters
        self.neighbors_dict = node_neighbors_dict
        self.triangles_list = triangles_list
        self.nodes = list(node_neighbors_dict.keys())
        self.N = len(node_neighbors_dict.keys())
        self.I = int(I_percentage * self.N/100)
        
        #Initial setup
        #I save the infected nodes of the first initialisation in case I want to repeat several runs with 
        #the same configuration
        self.initial_infected_nodes = self.initial_setup()
        
    def initial_setup(self, fixed_nodes_to_infect=None, print_status=True):
        #going to use this to store the agents in each state
        self.sAgentSet = set()
        self.iAgentSet = set()
        
        #and here we're going to store the counts of how many agents are in each
        #state @ each time step
        self.iList = []
        
        self.t = 0
        
        #start with everyone susceptible
        for n in self.nodes:
            self.sAgentSet.add(n)
        
        #infect nodes
        if fixed_nodes_to_infect==None: #the first time I create the model (the instance __init__)
            infected_this_setup=[]
            for ite in range(self.I): #we will infect I agents
                #select one to infect among the supsceptibles
                to_infect = random.choice(list(self.sAgentSet))
                self.infectAgent(to_infect)
                infected_this_setup.append(to_infect)
        else: #I already have run the model and this is not the first run, I want to infect the same nodes
            infected_this_setup=[]
            for to_infect in fixed_nodes_to_infect:
                self.infectAgent(to_infect)
                infected_this_setup.append(to_infect)
        #if print_status: print 'Setup:', self.N, 'nodes', self.I, 'infected'
        return infected_this_setup
    
    def infectAgent(self,agent):
        self.iAgentSet.add(agent)
        self.sAgentSet.remove(agent)
        return 1
    
    def recoverAgent(self,agent):
        self.sAgentSet.add(agent)
        self.iAgentSet.remove(agent)
        return -1
    
    def run(self, t_max, beta1, beta2, mu, print_status):
        self.t_max = t_max
        
        while len(self.iAgentSet) > 0 and len(self.sAgentSet) != 0 and self.t<=self.t_max:
            newIlist = set()
            
            #STANDARD CONTAGION
            #we only need to loop over the agents who are currently infectious
            for iAgent in self.iAgentSet:
                #expose their network neighbors
                for agent in self.neighbors_dict[iAgent]: 
                    #given that the neighbor is susceptible
                    if agent in self.sAgentSet:
                        #infect it with probability beta1
                        if (random.random() <= beta1): 
                            newIlist.add(agent)
              
            #TRIANGLE CONTAGION
            for triangle in self.triangles_list:
                n1, n2, n3 = triangle
                if n1 in self.iAgentSet:
                    if n2 in self.iAgentSet:
                        if n3 in self.sAgentSet:
                            #infect n3 with probability beta2
                            if (random.random() <= beta2): 
                                newIlist.add(n3)
                    else:
                        if n3 in self.iAgentSet:
                            #infect n2 with probability beta2
                            if (random.random() <= beta2): 
                                newIlist.add(n2)
                else:
                    if (n2 in self.iAgentSet) and (n3 in self.iAgentSet):
                        #infect n1 with probability beta2
                        if (random.random() <= beta2): 
                            newIlist.add(n1)
                
            #Update only now the nodes that have been infected
            for n_to_infect in newIlist:
                self.infectAgent(n_to_infect)
            
            
            #for recoveries
            newRlist = set()
            #In case all the individuals are infected I have to stop the recovery process. So I do the
            #recovery only if there is at least one individual not infected
            if len(self.iAgentSet)<self.N:
            
                for recoverAgent in self.iAgentSet:
                    #if the agent has just been infected it will not recover this time
                    if recoverAgent in newIlist:
                        continue
                    else:
                        if (random.random() <= mu): 
                            newRlist.add(recoverAgent)

            #Update only now the nodes that have been infected
            for n_to_recover in newRlist:
                self.recoverAgent(n_to_recover)
            
            #then track the number of individuals in each state
            self.iList.append(len(self.iAgentSet))
            
            #increment the time
            self.t += 1

        #and when we're done, return all of the relevant information
        if print_status: print('beta1', beta1, 'Done!', len(self.iAgentSet), 'infected agents left')

        return self.iList
    
    def get_stationary_rho(self, normed=True, last_k_values = 100):
        i = self.iList
        if len(i)==0:
            return 0
        if normed:
            i = 1.*np.array(i)/self.N
        if i[-1]==1:
            return 1
        elif i[-1]==0:
            return 0
        else:
            avg_i = np.mean(i[-last_k_values:])
            avg_i = np.nan_to_num(avg_i) #if there are no infected left nan->0   
            return avg_i

def run_one_simulation(args):
    
    it_num, N, p1, p2, lambda1s, lambdaD_target, I_percentage, t_max, mu  = args
    print('It %i initialized'%it_num)
    
    node_neighbors_dict, triangles_list = generate_my_simplicial_complex_d2(N,p1,p2)

    real_k = 1.*sum([len(v) for v  in node_neighbors_dict.values()])/len(node_neighbors_dict)
    real_kD = 3.*len(triangles_list)/len(node_neighbors_dict)

    print('It %i, created SC with k1=%.1f and k2=%.1f'%(it_num,real_k,real_kD))
    
    beta1s = 1.*(mu/real_k)*lambda1s
    beta2 = 1.*(mu/real_kD)*lambdaD_target
    
    rhos = [] #here I'll store the rho(t)
    
    for beta1 in beta1s:
        mySimplagionModel = SimplagionModel(node_neighbors_dict, triangles_list, I_percentage)
        mySimplagionModel.initial_setup(fixed_nodes_to_infect = mySimplagionModel.initial_infected_nodes);
        results = mySimplagionModel.run(t_max, beta1, beta2, mu, print_status=False)
        rho = mySimplagionModel.get_stationary_rho(normed=True, last_k_values = 100)
        rhos.append(rho)

    print('It %i, simulation has finished'%(it_num))

    return rhos, real_k, real_kD

    
def generate_my_simplicial_complex_d2(N,p1,p2):
    
    """Our model"""
    
    #I first generate a standard ER graph with edges connected with probability p1
    G = nx.fast_gnp_random_graph(N, p1, seed=None)
    
    if not nx.is_connected(G):
        giant = list(nx.connected_components(G))[0]
        G = nx.subgraph(G, giant)
        print('not connected, but GC has order %i ans size %i'%(len(giant), G.size())) 

    triangles_list = []
    G_copy = G.copy()
    
    #Now I run over all the possible combinations of three elements:
    for tri in combinations(list(G.nodes()),3):
        #And I create the triangle with probability p2
        if random.random() <= p2:
            #I close the triangle.
            triangles_list.append(tri)
            
            #Now I also need to add the new links to the graph created by the triangle
            G_copy.add_edge(tri[0], tri[1])
            G_copy.add_edge(tri[1], tri[2])
            G_copy.add_edge(tri[0], tri[2])
            
    G = G_copy
             
    #Creating a dictionary of neighbors
    node_neighbors_dict = {}
    for n in list(G.nodes()):
        node_neighbors_dict[n] = G[n].keys()           
                
    #print len(triangles_list), 'triangles created. Size now is', G.size()
        
    #avg_n_triangles = 3.*len(triangles_list)/G.order()
    
    #return node_neighbors_dict, node_triangles_dict, avg_n_triangles
    #return node_neighbors_dict, triangles_list, avg_n_triangles
    return node_neighbors_dict, triangles_list

def get_p1_and_p2(k1,k2,N):
    p2 = (2.*k2)/((N-1.)*(N-2.))
    p1 = (k1 - 2.*k2)/((N-1.)- 2.*k2)
    if (p1>=0) and (p2>=0):
        return p1, p2
    else:
        raise ValueError('Negative probability!')
        
def find_cut(rhos_array):
    #First index with non-zero value >1
    cut = min(np.argwhere(np.count_nonzero(rhos_array, axis=0)>1))[0]
    return cut

def parse_results(results, cut):
    
    rhos_array, real_k_list, real_kD_list = [], [], []
    
    for rhos, real_k, real_kD in results:
        real_k_list.append(real_k)
        real_kD_list.append(real_kD)
        rhos_array.append(rhos)
        
    rhos_array = np.array(rhos_array)
    real_kD_list = np.array(real_kD_list)    
    real_k_list = np.array(real_k_list)    
    
    avg_kD = real_kD_list.mean(axis=0)
    avg_k = real_k_list.mean(axis=0)

    if cut==False:    
        avg_rhos = np.mean(rhos_array, axis=0)
        avg_kD = real_kD_list.mean(axis=0)
        avg_k = real_k_list.mean(axis=0)

        #std_rhos = np.std(rhos_array, axis=0)
        return avg_rhos, avg_k, avg_kD    
    
    else:
        cut_point = find_cut(rhos_array)
        cut_rhos_array = []
        
        for rhos, _, _ in results:
            clean_rhos = []
            for i, rr in enumerate(rhos):
                if i<cut_point:
                    clean_rhos.append(rr)
                elif rr==0:
                    clean_rhos.append(np.nan)
                else:
                    clean_rhos.append(rr)
            cut_rhos_array.append(clean_rhos) 
        
        cut_rhos_array = np.array(cut_rhos_array)
        avg_rhos = np.nanmean(cut_rhos_array, axis=0)
        #std_rhos = np.nanstd(rhos_array, axis=0)
        
        return avg_rhos, avg_k, avg_kD    
    
#Function for MF
def get_rho_MF(l, lD):
    rho1 = (lD-l + np.sqrt((l-lD)**2 - 4.*lD*(1-l)))/(2*lD)
    if rho1>0: 
        return rho1
    else:
        return 0
