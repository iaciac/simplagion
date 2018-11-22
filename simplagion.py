
# coding: utf-8

"""A Python function that creates a random simplicial complex,
    as defined in Iacopini, I., Petri, G., Barrat, A., & Latora, V. (2018).
    Simplicial models of social contagion. arXiv preprint arXiv:1810.07031.
"""

#    Copyright (C) 2018 by
#    Iacopo Iacopini <i.iacopini@qmul.ac.uk>
#    All rights reserved.
#    BSD license.

import networkx as nx
import itertools
import random

__all__ = ['generate_my_simplicial_complex_d2']


def generate_my_simplicial_complex_d2(N,k1,k2):
    """Returns a random simplicial complex of dimension d=2 (see [1]).
        
        Args
        ----
        N: int
        Number of nodes (0-simplices).
        
        k1: int
        Expected average 1-degree (number of 1-simplices incident on a node)
        
        k2: int
        Expected average 2-degree (number of 2-simplices incident on a node)
        
        Returns
        -------
        node_neighbors_dict: dict
        Dictionary having nodes as keys and lists of neighbors as values.
        
        triangles_list: list
        List of 2-simplices generated [(i,j,k),...]
        
        avg_n_triangles: float
        Average number of 2-simplices
                
        References
        ----------
        .. [1] Iacopini, I., Petri, G., Barrat, A., & Latora, V. (2018).
        "Simplicial models of social contagion".
        arXiv preprint arXiv:1810.07031.
        
        """
        
    #Calculating - approximated - probabilities
    p2 = (2.*k2)/((N-1.)*(N-2.))
    p1 = (k1 - 2.*k2)/(N-1.)
    if (p1<0) or (p2<0): raise ValueError('Negative probability!')
    
    #I first generate a standard ER graph with edges connected with probability p1
    G = nx.fast_gnp_random_graph(N, p1, seed=None)
    
    triangles_list = []
    
    #Now I run over all the possible combinations of three elements:
    for tri in itertools.combinations(list(G.nodes()),3):
        #And I create the triangle with probability p2
        if random.random() <= p2:
            #I close the triangle.
            triangles_list.append(tri)
            
            #Now I also need to add the new links to the graph created by the triangle
            G.add_edge(tri[0], tri[1])
            G.add_edge(tri[1], tri[2])
            G.add_edge(tri[0], tri[2])
     
    if not nx.is_connected(G):
        giant = max(nx.connected_component_subgraphs(G), key=len)
        print 'not connected, but GC has order ', giant.order(), 'and size', giant.size()
        G = giant
    
    #Creating a dictionary of neighbors
    node_neighbors_dict = {}
    for n in G.nodes():
        node_neighbors_dict[n] = G[n].keys()
    
    print len(triangles_list), 'triangles created. Size now is', G.size()
        
    avg_n_triangles = 3.*len(triangles_list)/G.order()
    
    return node_neighbors_dict, triangles_list, avg_n_triangles

