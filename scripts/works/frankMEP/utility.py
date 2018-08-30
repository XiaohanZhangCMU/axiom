""" A set of utility functions to process set operations of lists of nucleus indices
"""
import numpy as np
import pickle
#import matplotlib.pyplot as plt

def removeNucleusAtoms(cfg, nucleus, NP):
    idx = np.setdiff1d(np.xarrange(NP), nucleus)
    cfg = cfg[idx,:]
    return cfg

# Return a list of atom adjacent to, but not in, nucleus
def find_nbr_atoms(nbrlist, nucleus, totIdx):
    B = nbrlist[nucleus,:]
    A = np.setdiff1d(np.extract(B>=0,B), nucleus)
    return np.intersect1d(totIdx, A)

def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)

def merge_obj(obj_1, obj_2):
    obj = obj_1.copy()
    obj.update(obj_2)
    return obj

def nucleus2bits(nucleus, totIdx):
    bits = np.in1d(totIdx, nucleus)
    bits = bits.astype(int)
    return bits

def bits2nucleus(bits, totIdx):
    nucleus = totIdx[np.where(bits==1)]
    nucleus = nucleus.astype(int)
    return nucleus

def bits2str(bits):
    mybits=""
    for k in bits:
      mybits+=str(k)
    return mybits

def str2bits(bitstr):
    return np.array(map(int, bitstr))

def getnbrlist(nucleus, nbrlist):
    tmp = nbrlist[nucleus,:]
    return  np.extract(tmp>=0,tmp)

def nbitsdiff(bit1, bit2):
    #  print('diff by {0}'.format(len((np.where((bit1^bit2)==1))[0])))
    return len((np.where((bit1^bit2)==1))[0])

def nbits1(bit):
    return len(np.where(bit==1)[0])

