import networkx as nx
import string
import json
import numpy as np
from time import time

def extract_networks(data_dir, dataset, n_minutes=5,original_nets=True):
    f = open(data_dir+'/tij_' + dataset +'.dat')
    (t0,i,j) = map(int,string.split(f.readline()))
    if dataset not in ['LyonSchool','LH10']:
        t0 = t0*20
    f.close()
    delta_t = 20*3*n_minutes;   
    if original_nets==True:
        originalnetworks = {}
    aggnetworks = {}
    f = open(data_dir+'/tij_' + dataset +'.dat')
    for line in f:
        (t,i,j) = map(int,string.split(line))
        if dataset not in ['LyonSchool','LH10']:
            t = t*20
        if original_nets==True:
            if t not in originalnetworks:
                originalnetworks[t] = nx.Graph()
            originalnetworks[t].add_edge(i,j)
        #this is a trick using the integer division in python
        aggtime = t0 + ((t-t0)/delta_t)*delta_t 
        if aggtime not in aggnetworks:
            aggnetworks[aggtime] = nx.Graph()
        aggnetworks[aggtime].add_edge(i,j)
    f.close();
    if original_nets==True:
        return originalnetworks, aggnetworks;
    else:
        return aggnetworks;
    
def extract_cliques(gs):
    listsaggcliques = {}
    for t in sorted(gs.keys()):
        listsaggcliques[t] = list(nx.find_cliques(gs[t]));
    return listsaggcliques;
    
def clique_weights(cliques):
    from collections import Counter;
    tot_c = [];
    for t in cliques:
        tot_c.extend(map(frozenset,cliques[t]))
    return Counter(tot_c);

def average_clique_size(ws):
    return np.sum(map(lambda x: 1.0 * ws[x] * len(x), ws.keys()))/np.sum(ws.values());

def clean_non_maximal(ws):
    sd = dict(zip(ws.keys(), map(len,ws.keys())));
    import operator
    sizes = set(map(len,ws.keys()));
    sorted_sd = sorted(sd.items(), key=operator.itemgetter(1));
    simplices = dict.fromkeys(list(sizes),[]);
    maximal_simplices = {};
    for x in ws:
        maximal = True;
        for xx in ws:
            if (len(x)<len(xx)):
                if (set(x)<set(xx)):
                    maximal=False;
                    break;
        if maximal:
            maximal_simplices[x] = ws[x];
    return maximal_simplices;

def save_cliques(ws, data_dir, dataset,n_minutes, thr=None):
    if thr==None:
        ls = map(list,ws.keys());
    else:
        ls = [list(x) for x in ws if ws[x]>=thr];
    jd = open(data_dir+'aggr_'+str(n_minutes)+'min_cliques_'+dataset+'.json','w')
    json.dump(ls,jd)
    jd.close()
    return;


def reweighting_all_cliques(mc):
    from itertools import combinations;
    cd = [];
    for c in mc:
        for d in range(2, len(c)+1):
            cd.extend(list(combinations(c,d)));
    cd = list(set(cd));
    cd = dict.fromkeys(cd,0);
    for c in mc:
        for d in range(2, len(c)+1): #I weigh only from the edges upwards
            for comb in combinations(c,d):
                cd[comb]+=1;
    return cd;
            
def quantile_cut(rwm,quant):
    q = np.quantile(rwm.values(), quant);
    return [k for k in rwm if rwm[k]>=q];
    
def limit_dimension(rwm,d):
    return [k for k in rwm if len(k)<=d];


