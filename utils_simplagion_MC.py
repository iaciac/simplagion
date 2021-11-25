import numpy as np
import json
import random
import networkx as nx
from itertools import combinations
from collections import defaultdict


def find_cut(results):
    rhos_array = np.array(results)
    #First index with non-zero value >1
    cut = min(np.argwhere(np.count_nonzero(rhos_array, axis=0)>1))[0]
    return cut

def parse_results(results, cut):
    
    if cut==False:
        rhos_array = np.array(results)
        avg_rhos = np.mean(rhos_array, axis=0)
        #std_rhos = np.std(rhos_array, axis=0)
        return avg_rhos#, std_rhos
    else:
        cut_point = find_cut(results)

        rhos_array = []
        for rhos in results:
            clean_rhos = []
            for i, rr in enumerate(rhos):
                if i<cut_point:
                    clean_rhos.append(rr)
                elif rr==0:#rr<0.1:#rr==0:
                    clean_rhos.append(np.nan)
                else:
                    clean_rhos.append(rr)
            rhos_array.append(clean_rhos) 
        
        rhos_array = np.array(rhos_array)
        avg_rhos = np.nanmean(rhos_array, axis=0)
        std_rhos = np.nanstd(rhos_array, axis=0)
        
        return np.nan_to_num(avg_rhos), np.nan_to_num(std_rhos)
    
def get_tri_neighbors_dict(triangles_list):
    tri_neighbors_dict = defaultdict(list)
    for i, j, k in triangles_list:
        tri_neighbors_dict[i].append((j,k))
        tri_neighbors_dict[j].append((i,k))
        tri_neighbors_dict[k].append((i,j))
    return tri_neighbors_dict

def import_sociopattern_simcomp_SCM(dataset_dir, dataset, n_minutes, thr):
    filename = dataset_dir+'random_'+str(n_minutes)+'_'+str(thr)+'min_cliques_'+dataset+'.json'
    SCM_cliques_list = json.load(open(filename,'r'))
    
    #considering one realization of the SCM
    realization_number = random.choice(range(len(SCM_cliques_list)))
    cliques = SCM_cliques_list[realization_number] 
    node_neighbors_dict, triangles_list = create_simplicial_complex_from_cliques(cliques)
    
    N = len(node_neighbors_dict.keys())
    avg_k1 = 1.*sum([len(v) for v in node_neighbors_dict.values()])/N
    avg_k2 = 3.*len(triangles_list)/N 
    #ass = nx.degree_assortativity_coefficient(facet_list_to_graph(cliques))

    #print N, avg_k1, avg_k2, "Assortativity %.2f"%ass
    
    return node_neighbors_dict, triangles_list, avg_k1, avg_k2

def create_simplicial_complex_from_cliques(cliques):
    
    G = nx.Graph()
    triangles_list = set() #will contain list of triangles (2-simplices)
    
    for c in cliques:
        d = len(c)
        
        if d==2:
            i, j = c
            G.add_edge(i, j)
        
        elif d==3:
            #adding the triangle as a sorted tuple (so that we don't get both ijk and jik for example)
            triangles_list.add(tuple(sorted(c)))
            #adding the single links
            for i, j in combinations(c, 2):
                G.add_edge(i, j)
            
        else: #d>3, but I only consider up to dimension 3
            #adding the triangles
            for i, j, k in combinations(c, 3):
                triangles_list.add(tuple(sorted([i,j,k])))

            #adding the single links
            for i, j in combinations(c, 2):
                G.add_edge(i, j)
                
    if nx.is_connected(G)==False:
        print('not connected')
                
    #Creating a dictionary of neighbors
    node_neighbors_dict = {}
    for n in G.nodes():
        node_neighbors_dict[n] = G[n].keys()
        
    #converting the triangle set of tuples into a triangle list of lists
    triangles_list = [list(tri) for tri in triangles_list]

    return node_neighbors_dict, triangles_list

