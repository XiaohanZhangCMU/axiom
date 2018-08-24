import plotly.plotly as py
from plotly.graph_objs import *
from plotly.offline import plot
import networkx as nx
from visualize_network import visualize_html
from utility import writedata,writecncfg,save_obj,load_obj,bits2nucleus,str2bits,bits2str,nucleus2bits,removeNucleusAtoms,getnbrlist,nbitsdiff
from setget_db import select_atoms_within_ecllipse, select_totalatoms_on_slice, readnn, get_nbrlist, find_bdy_atoms, make_pairs, calculate_frk_energy, calculate_frk_energy, load_status_0, setget_db
import numpy as np

'''
Find minimum energy. Assign ID to each key
'''
dirname = '/home/xzhang11/Planet/Libs/MD++.git/runs/frankMEP/test_db_0.05/'
#dirname = 'runs/frankMEP/test_shortest_path_db_0.05/'
db = load_obj(dirname+'db_0.05')
#db = load_obj(dirname+'test_single_graph_db')

ID = { }
db_simple = { } 
id = 0;
min_energy = 1e5;
for key, val in db.iteritems():
    ID[key] = id
    id = id + 1
    if val < min_energy:
        min_energy = val
    db_simple[key] = val

print("min_energy = {0}".format(min_energy))
#Find configuration of start (fewest atoms) and end (most atoms)

'''
Assign start end points for search. start: fewest atoms. end: most
'''
idfewest = 0
idmost = 0
fewest = 1e5
most = -1e5
for key, val in db.iteritems():
    sz = len(np.where(str2bits(key)==1)[0]) #number of `1' in the bit array
    if sz < fewest:
        fewest = sz
	idfewest = ID[key]
#    else: 
    if sz > most:
        most = sz
        idmost = ID[key]

fews = []
mosts = []
energy_fews = []
energy_mosts = []
for key, val in db.iteritems():
    sz = len(np.where(str2bits(key)==1)[0])
    if sz == fewest:
        fews.append(ID[key])
	energy_fews.append(db[key])
    if sz == most:
        mosts.append(ID[key])
        energy_mosts.append(db[key])
	
idmost = mosts[energy_mosts.index(min(energy_mosts)) ]
idfewest = fews[energy_fews.index(min(energy_fews)) ]
print('idmost = {0}'.format(idmost))
print('idfewest = {0}'.format(idfewest))
print('most = {0}'.format(most))
print('fewest = {0}'.format(fewest))
print('most = {0}'.format(mosts))
print('fews = {0}'.format(fews))

'''
Read in graph and look for shortest path

'''
G = nx.read_gpickle(dirname+"/nx_graph.gpickle")
path =nx.shortest_path(G, source=idfewest, target=idmost) 
print(path)

nucleus = []
potential = []
with open(dirname+'shortest_path_potential.dat', 'w') as fp:
    id0 = 0
    (totIdx, atomIdx,cfg,H,normal,d,pt0,nbrlist,slice_nbrlist_u,slice_nbrlist_d, pairs)= load_status_0(dirname) 
    for path_id in path:
        for key, val in ID.iteritems():
            if val == path_id:
                nucleus.append(str2bits(key))
    	        potential.append(db[key])
    	        fp.write(str(db[key])+"\n")
                writecncfg(cfg[bits2nucleus(str2bits(key), totIdx),:], H, dirname+"path-"+str(id0))
		id0 += 1
	

'''
Find MST of G, (supposed to be MEP?)
'''
T = nx.minimum_spanning_tree(G)
P = sorted(T.edges(data=True))

print(P)


#quit()


'''
Find all simple path starting from the shortest one
'''

simple_paths = nx.shortest_simple_paths(G, source=idfewest, target=idmost)
#sz = 0
#for path in simple_paths:
#    print('sz = {0}'.format(sz))
#    sz += 1
print("size of simple_paths = {0}".format(sz))
with open(dirname+'all_simple_path_potential.dat', 'w') as fp:
    (totIdx, atomIdx,cfg,H,normal,d,pt0,nbrlist,slice_nbrlist_u,slice_nbrlist_d, pairs)= load_status_0(dirname) 
    id_path = 0
    for path in simple_paths:
        nucleus = []
        potential = []
        id0 = 0
	if len(path) > 28:
	    continue
        for path_id in path:
            for key, val in ID.iteritems():
                if val == path_id:
                    nucleus.append(str2bits(key))
        	    potential.append(db[key])
        	    fp.write(str(db[key])+" ")
                    writecncfg(cfg[bits2nucleus(str2bits(key), totIdx),:], H, dirname+"path-"+str(id_path)+"-"+str(id0))
    		    id0 += 1
	fp.write("\n")
	fp.flush()
	id_path+=1
    
#print(nucleus)
print(potential)
    

