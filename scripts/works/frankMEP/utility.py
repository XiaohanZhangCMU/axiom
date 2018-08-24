import os, sys
import numpy as np
import random
import pickle
from subprocess import call
from numpy import linalg as LA
from bitarray import bitarray
import matplotlib.pyplot as plt

# Return atom Id within an ecllipse given _SR
# pt1: center coordinate in _SR coordinate
# a, b: width and length
# normal: of the plane ecllipse lies on

def writedata(idx, data, finalcnfile):
  with open(finalcnfile+".dat", 'w') as the_file:
    for i in range(data.shape[0]):
      the_file.write(str(idx[i]) + " ")
      for j in range(3):
        the_file.write(str(data[i,j]) + " ")
      the_file.write('\n') 

def writecncfg(data, H, finalcnfile):
  with open(finalcnfile+".cn", 'w') as the_file:
    the_file.write(str(data.shape[0])+'\n')
    for i in range(data.shape[0]):
      for j in range(3):
        the_file.write(str(data[i,j]) + "\t")
      the_file.write('\n') 
    for i in range(3):
      for j in range(3):
        the_file.write(str(H[i,j]) + "\t")
      the_file.write('\n') 
    the_file.write('1 Si\n 0 0')
  call([ "sw_mc2_mpich",  "scripts/work/frankMEP/disl_nuc_hetero.tcl", "100", finalcnfile ])
  
def removeNucleusAtoms(cfg, nucleus, atomIdx):
  idx = np.setdiff1d(atomIdx, nucleus)
  cfg = cfg[idx,:]
  return cfg
  
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

def visualize_path(path):
  (X,Y,Z) = order_path(path)
  plt.plot(X,Y,'*-')
  plt.draw()

def write_path(path, name):
  (N,Y,X) = order_path(path)
  with open(name, "w") as fp:
    for val in N:
      fp.write(str(val)+" ")
    fp.write("\n")
    for val in Y:
      fp.write(str(val)+" ")
    fp.write("\n")
  
def order_path(path):
  X = np.array(list(path.keys()))
  Y = np.array(list(path.values()))
  N = [ ]
  for item in X:
    N.append(nbits1(str2bits(item)))
  N = np.array(N) 
  inds = N.argsort()
  Y = Y[inds]
  N = N[inds]
  X = X[inds]
  return (N, Y, X) 



