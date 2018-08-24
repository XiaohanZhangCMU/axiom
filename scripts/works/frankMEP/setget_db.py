import os, sys
import numpy as np
import pickle
from subprocess import call
from numpy import linalg as LA
from bitarray import bitarray
from utility import writedata,writecncfg,save_obj,load_obj,bits2nucleus,bits2str,nucleus2bits,removeNucleusAtoms,getnbrlist

def select_atoms_within_ecllipse(pt0, a, b, dirname):
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, slice_nbrlist_u, slice_nbrlist_d,pairs)= load_status_0(dirname) 

  print("Select nucleus within an elipse.")
  e0 = [1,0,0]; e1 = [0,1,0];                e2 = [0,0,1];
  ep0 = e1;     ep1 = np.cross(normal, ep0); ep2 = normal; 
  # coordinate transform
  Q = np.zeros((3,3))
  Q[0][0] = np.dot(ep0,e0)
  Q[1][1] = np.dot(ep1,e1)
  Q[2][2] = np.dot(ep2,e2)
  Q[0][1] = np.dot(ep0, e1)
  Q[0][2] = np.dot(ep0, e2)
  Q[1][0] = np.dot(ep1, e0)
  Q[1][2] = np.dot(ep1, e2)
  Q[2][0] = np.dot(ep2, e0)
  Q[2][1] = np.dot(ep2, e1)

  nucleus = []; 
  for i in range(cfg.shape[0]):
    _x0 = cfg[i,:] - pt0; # origin shifts to the picked point
    _x = np.dot(Q,_x0);
    x = _x[0]; y = _x[1]; z = _x[2]; 
    if (x)*(x)/a/a + (y)*(y)/b/b <= 1 :
      nucleus.append(i);
  nucleus = np.array(nucleus)
  return np.intersect1d(nucleus, totIdx)
  
# Find atom id that is eps away from a plane defined by normal and _pt
# If opt > 0, choose atoms above the plane, otherwise below.
def select_totalatoms_on_slice(eps, _pt, cfg, normal, atomIdx, opt):
  x0 = _pt[0]; y0 = _pt[1]; z0 = _pt[2];
  d = normal[0] * x0 + normal[1] * y0 + normal[2] * z0;
  D = np.ones(cfg[:,0].shape) * d ; 
  cond = np.abs(normal[0] * cfg[:,0] + normal[1] * cfg[:,1]  + normal[2] * cfg[:,2] -D) < eps ;
  idx = np.extract(cond,atomIdx)
  if opt != 0:
    D = np.ones(cfg[idx,0].shape) * d ; 
    if opt >0:
      cond = normal[0] * cfg[idx,0] + normal[1] * cfg[idx,1]  + normal[2] * cfg[idx,2] -D > 0 ;
    elif opt <0:
      cond = normal[0] * cfg[idx,0] + normal[1] * cfg[idx,1]  + normal[2] * cfg[idx,2] -D < 0 ;
    idx = np.extract(cond, idx)
  return idx


# Build nbrlist for a group of atoms (totIdx)
def readnn(cfg, maxnbrs, nbrlistname):
  print("Read nnlist dumped by MD++ status = 0 tcl")
  nbrlist = -1* np.ones((cfg.shape[0], maxnbrs)) 
  with open(nbrlistname, 'r') as f:
    lines = f.readlines()
    I = 0
    for line in lines:
      array = np.fromstring(line, dtype=int, sep = ' ')
      J = 0
      for id in array:
        nbrlist[I, J] = id
        J += 1
      I += 1

  return nbrlist

# Build nbrlist for a group of atoms (totIdx)
def get_nbrlist(totIdx, cfg, maxnbrs, cutoff, nbrlistname, opt):
  if not opt:
    nbrlist = np.load(nbrlistname)
  else:
    print("Buidling neighbor list for atoms on slice.")
    print("Patience......")
    nbrlist = -1* np.ones((cfg.shape[0], maxnbrs)) 
    nbrlist = nbrlist.astype(int)
    for atom_I in totIdx:
      cnt = 0
      for atom_J in totIdx:
        if atom_I != atom_J:
          a = LA.norm(cfg[atom_I, :]-cfg[atom_J,:], 2) 
          if a < cutoff:
            nbrlist[atom_I, cnt] = atom_J
            cnt += 1
            assert (cnt < maxnbrs), "Need to increase maxnbrs."
            assert (cnt > 0), "Need to increase cutoff."
    np.save(nbrlistname, nbrlist)
  return nbrlist

# Return a list of atom ID that is next but not within the nucleus
def find_bdy_atoms(nbrlist, nucleus, totIdx, cfg):
  B = nbrlist[nucleus,:]
  #print("B.size = {0}".format(B.size))
  A = np.setdiff1d(np.extract(B>=0,B), nucleus)
  #print("A.size = {0}".format(A.size))
  #print(A)
  #(totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, slice_nbrlist_u, slice_nbrlist_d, pairs)= load_status_0('runs/frankMEP/test_db_0.05/') 
  #writecncfg(cfg[A,:], H, "tmp")
  #print("totIdx.size = {0}".format(totIdx.size))
  #print(totIdx)
  #print(np.intersect1d(totIdx, A))
  return np.intersect1d(totIdx, A)

def make_pairs(totIdx_1, totIdx_2, cfg):
  if len(totIdx_1) < len(totIdx_2):
    t1 = totIdx_1
    t2 = totIdx_2
  else:
    t1 = totIdx_2
    t2 = totIdx_1
  pairs = [ ]

  for atom_I in t1:
    mindis = 1e5
    for atom_J in t2:
      tmp = LA.norm(cfg[atom_I,:]-cfg[atom_J,:]) 
      if tmp < mindis:
        mindis = tmp
	minatom = atom_J
    pairs.append([atom_I, minatom])
  return np.array(pairs).astype(int)

# Call tcl file to remove nucleus from 0K_0.0_relaxed_001.cn
# Then relax and compute energy which is written to EPOT_2.dat
def calculate_frk_energy(nucleus, strain, dirname, id0):
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, slice_nbrlist_u, slice_nbrlist_d, pairs)= load_status_0(dirname) 
  #print nucleus
  writedata(nucleus, cfg[nucleus,:],  dirname+"nucleus-"+str(id0))

  idx_u = np.intersect1d(getnbrlist(nucleus, nbrlist), slice_nbrlist_u)
  idx_d = np.intersect1d(getnbrlist(nucleus, nbrlist), slice_nbrlist_d)
  writedata(idx_u, cfg[idx_u,:],  dirname+"nucleus_nbrlist_u-"+str(id0))
  writedata(idx_d, cfg[idx_d,:],  dirname+"nucleus_nbrlist_d-"+str(id0))
  writecncfg(cfg[idx_u,:], H,  dirname+"nucleus_nbrlist_u-"+str(id0))
  writecncfg(cfg[idx_d,:], H,  dirname+"nucleus_nbrlist_d-"+str(id0))

  writecncfg(cfg[nucleus,:], H,  dirname+"nucleus-"+str(id0))
  call([ "sw_mc2_mpich","scripts/work/frankMEP/disl_nuc_hetero.tcl","11",str(id0),str(normal[0]), str(normal[1]), str(normal[2]), str(d), str(pt0[0]), str(pt0[1]), str(pt0[2]), dirname, strain]);
  while not os.path.exists(dirname+'EPOT_2.dat'):
    continue
  energy = np.loadtxt(dirname+'EPOT_2.dat')
  return energy

# Load all information from status == 0, including:
# System: cfg, H, atomIdx (simply 1:NP, may not need it)
# Slice: A*x+B*y+C*z-D=0
# Atom id of the whole slice: totIdx
def load_status_0(dirname):
  totIdx = np.load(dirname+"totIdx.npy")
  atomIdx = np.load(dirname+"atomIdx.npy")
  cfg = np.load(dirname+"cfg.npy")
  H = np.load(dirname+"H.npy")
  ctrlparams = np.load(dirname+"ctrlparams.npy")
  normal= ctrlparams[0:3]; 
  d = ctrlparams[3];
  pt0 = ctrlparams[4:7]
  nbrlist = np.load(dirname+"nbrlist.npy")
  slice_nbrlist_u = np.load(dirname+"slice_nbrlist_u.npy")
  slice_nbrlist_d = np.load(dirname+"slice_nbrlist_d.npy")
  pairs = np.load(dirname+"pairs.npy")
  return (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, slice_nbrlist_u, slice_nbrlist_d, pairs)

def setget_db(nucleus, strain, cohesv, dirname, id0):
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, slice_nbrlist_u, slice_nbrlist_d,pairs)= load_status_0(dirname) 

  db_file = dirname+"db_"+strain+".2rm";
  if os.path.exists(db_file+".pkl"):
    db  = load_obj(db_file)
  else:
    db = { }

  bitarr = nucleus2bits(nucleus, totIdx)
  bitstr = bits2str(bitarr)
  if bitstr in db:
    print("data base has = {0} data points".format(len(db)))
    return db[bitstr]
  else:
    energy = calculate_frk_energy(nucleus, strain, dirname, id0)
    energy -= cohesv * nucleus.size;
    db[bitstr] = energy
    save_obj(db, db_file)
    return energy

def test_single_graph_setget_db(nucleus, energy, dirname):
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, slice_nbrlist_u, slice_nbrlist_d,pairs)= load_status_0(dirname) 

  db_file = dirname+"test_single_graph_db";
  if os.path.exists(db_file+".pkl"):
    db  = load_obj(db_file)
  else:
    db = { }

  bitarr = nucleus2bits(nucleus, totIdx)
  bitstr = bits2str(bitarr)
  db[bitstr] = energy
  save_obj(db, db_file)
  return energy

