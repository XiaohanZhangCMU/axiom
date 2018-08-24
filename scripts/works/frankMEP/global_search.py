import os, sys
import numpy as np
import pickle
import math
import random
from subprocess import call
from numpy import linalg as LA
from bitarray import bitarray
from utility import writedata,writecncfg,save_obj,load_obj,bits2nucleus,bits2str, str2bits,nbits1, nucleus2bits,removeNucleusAtoms, getnbrlist, order_path, write_path
from setget_db import select_atoms_within_ecllipse, select_totalatoms_on_slice, readnn, get_nbrlist, find_bdy_atoms, make_pairs, calculate_frk_energy, calculate_frk_energy, load_status_0, setget_db

#Prepare initial conditions, dumping a number of npy files for other status
#In principle, you don`t need to rerun status = 1 if you have done it before
def status0(strain, dirname, cohesv): 
  #call([ "bin1/sw_mc2_mpich",  "scripts/work/frankMEP/disl_nuc_hetero.tcl", "0", "n", "001", "1", dirname, strain ])

  # Run status == 0 to pick the ctrl pts to define slice
  # atoms between these two planes will be removed.
  # pt1_up/pt1_d----pt2_up/pt2_d  --> left glide plane
  # pt3_up/pt3_d----pt4_up/pt4_d  --> right glide plane
  # pt0: middle of pt1 and pt3, thus on shuffle plane
  pt1d = np.array([0.0208,0, 0.3784]); pt1u = np.array([0.0208,0, 0.3911])
  pt2d = np.array([0.1042, 0, 0.2734]);pt2u = np.array([0.1042, 0, 0.2865])
  pt3d = np.array([0.0625,0, 0.3773]); pt3u = np.array([0.0625,0, 0.3911])
  pt4d = np.array([0.1250, 0, 0.2995]);pt4u = np.array([0.1250, 0, 0.3125])
  pt1 = 0.5*(pt1d + pt1u);  pt2 = 0.5*(pt2d + pt2u)
  pt3 = 0.5*(pt3d + pt3u);  pt4 = 0.5*(pt4d + pt4u)
  pt0 = 0.5*(pt1+pt3); 

  v = pt2- pt1
  v = v/LA.norm(v)
  normal = np.cross(v, [0,1,0])

  x0 = pt0[0]; y0 = pt0[1]; z0 = pt0[2];
  d = normal[0] * x0 + normal[1] * y0 + normal[2] * z0;
  ctrlparams =np.array( [  normal[0], normal[1], normal[2], d, pt0[0], pt0[1], pt0[2] ] )
  pt1pt3 = np.array([ pt1[0], pt1[1], pt1[2], pt3[0], pt3[1], pt3[2]])

  print("ctrlparams = {0}".format(ctrlparams))

  # Read in perfect configuration .cn file
  data=np.genfromtxt(dirname+"0K_0.0_relaxed_surf001.cn",skip_header=1,skip_footer=2);
  cfg  = data[:-3,:]
  H    = data[-3::,:]
  atomIdx = np.arange(cfg.shape[0])

  eps = 0.02
  cutoff = 0.07
  maxnbrs = 40

  assert(np.abs(np.sqrt(v[0]*v[0]+v[1]*v[1] + v[2]*v[2])- LA.norm(v))<1e-10)
  totIdx_1 = select_totalatoms_on_slice(eps, pt1, cfg, normal, atomIdx, 1)
  totIdx_2 = select_totalatoms_on_slice(eps, pt3, cfg, normal, atomIdx, -1)
  print("pt3 = {0}".format(pt3))
  print("pt1 = {0}".format(pt1))
  totIdx = np.union1d(totIdx_1, totIdx_2)
  pairs = make_pairs(totIdx_1, totIdx_2, cfg)

  # Build neighborhood list for totIdx
  nbrlist = readnn(cfg,maxnbrs,dirname+"nn.dat");
  nbrlist = nbrlist.astype(int)
  totIdx_u = select_totalatoms_on_slice(2*eps, pt0, cfg, normal, atomIdx, 1)
  totIdx_d = select_totalatoms_on_slice(2*eps, pt0, cfg, normal, atomIdx, -1)
  slice_nbrlist_u = np.intersect1d(getnbrlist(totIdx_1, nbrlist), totIdx_u)
  slice_nbrlist_d = np.intersect1d(getnbrlist(totIdx_2, nbrlist), totIdx_d)
  #nbrlist = get_nbrlist(totIdx, cfg, maxnbrs, cutoff, dirname+"nbrlist.npy", True);
  #nbrlist1 = nbr4all[totIdx,:]

  eps = 0.02 # how thick you want to sandwich the slice
#  idx_1 = findNearByAtoms(eps,pt1, cfg, normal, atomIdx)
#  idx_2 = findNearByAtoms(eps,pt3, cfg, normal, atomIdx)
  idx_1 = select_totalatoms_on_slice(eps,pt1, cfg, normal, atomIdx, 0)
  idx_2 = select_totalatoms_on_slice(eps,pt3, cfg, normal, atomIdx, 0)
  idx = np.union1d(idx_1,idx_2)

  writedata(idx, cfg[idx,:], dirname+"nbratoms")
  writecncfg(cfg[idx,:], H, dirname+"nbratoms")


  writecncfg( cfg[totIdx,:], H, dirname+"slice")
  writedata(totIdx, cfg[totIdx,:], dirname+"slice")

  writecncfg( cfg[slice_nbrlist_u,:], H, dirname+"slice_nbrlist_u")
  writecncfg( cfg[slice_nbrlist_d,:], H, dirname+"slice_nbrlist_d")

  writecncfg( cfg[pairs[:,0],:], H, dirname+"pairs_u")
  writecncfg( cfg[pairs[:,1],:], H, dirname+"pairs_d")


#  for i in range(pairs.shape[0]):
#    writecncfg( cfg[pairs[i,:],:], H, dirname+"pairs_"+str(i))

  writecncfg( cfg[totIdx_u,:], H, dirname+"totIdx_u")
  writecncfg( cfg[totIdx_d,:], H, dirname+"totIdx_d")

  writecncfg( cfg[totIdx_1,:], H, dirname+"totIdx_1")
  writecncfg( cfg[totIdx_2,:], H, dirname+"totIdx_2")

  np.save(dirname+"totIdx.npy", totIdx)
  np.save(dirname+"atomIdx.npy", atomIdx)
  np.save(dirname+"cfg.npy", cfg)
  np.save(dirname+"H.npy", H)
  np.save(dirname+"ctrlparams.npy", ctrlparams)
  np.save(dirname+"pt1pt3.npy", pt1pt3)
  np.save(dirname+"nbrlist.npy", nbrlist)
  np.save(dirname+"slice_nbrlist_u.npy", slice_nbrlist_u)
  np.save(dirname+"slice_nbrlist_d.npy", slice_nbrlist_d)
  np.save(dirname+"pairs.npy", pairs)
  exit(0)
 
#Populate data base with ellipse of different asp ratio.
#For each nucleus, tcl file is called to remove the atoms.
#Energy of the collapsed frank partial is then calculated.
def status1(strain, dirname, cohesv): 
  Nmax = 20

  pt1pt3 = np.load(dirname+"pt1pt3.npy")
  pt1 = pt1pt3[0:3]

  asps = np.linspace(0,1,25)
  # Some controls on ellipse loops
  step0 = 0 # An file identifier to tell configuration files apart
  for asp in np.linspace(0.1, 3, 30):
    a0 = 0.02; b0 = a0*asp; max_radius = 0.75

    for step in range(Nmax):
      # Specify ellipse loops
      nucleus = select_atoms_within_ecllipse(pt1, a0, b0, dirname)   
      if nucleus.size == 0:
        a0 += (max_radius-a0)/Nmax; b0 = a0*asp;
        continue

      # Find energy for frank partial with nucleus removed
      energy = setget_db(nucleus, strain, cohesv, dirname, step0)
      print("Energy of step {0} = {1}. a0 = {2}".format(step, energy, a0))
      with open(dirname+"potential.dat", "a") as fp:
        fp.write(str(energy)+"\n")

      a0 += (max_radius-a0)/Nmax; b0 = a0*asp;
      step0 += 1

#Read in a data base and convert each entry to write out MD configuration
def status2(strain, dirname, cohesv): 
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, slice_nbrlist_u, slice_nbrlist_d, pairs)= load_status_0(dirname) 
  db_file = dirname+"db_"+strain;
  if os.path.exists(db_file+".pkl"):
    db  = load_obj(db_file)
    print("data base has = {0} data points".format(len(db)))
  else:
    print("data base {0} not found".format(db_file))
  
  global PY3
  if sys.version_info[0] < 3:
    PY3 = False
    print("Use Python2")
  else:
    PY3 = True
    print("Use Python3")

  step0 = -10000 # this is a file identifier
#  if PY3:
#    for key, value in db.items():
  for key, value in db.iteritems():
    nucleus = bits2nucleus(str2bits(key), totIdx)
    energy = calculate_frk_energy(nucleus, strain, dirname, step0)
    step0 += 1

#Get a bumpy energy barrier curve for gradually bigger loops 
#plot with 'plot_energy_barrier.m'
def status3(strain, dirname, cohesv): 
  Nmax = 20
  pt1pt3 = np.load(dirname+"pt1pt3.npy")
  pt1 = pt1pt3[0:3]
  a0 = 0.02; b0 = a0; max_radius = 0.75

  step0 = 0 # An file identifier to tell configuration files apart
  for step in range(Nmax):
    # Specify ellipse loops
    nucleus = select_atoms_within_ecllipse(pt1, a0, b0, dirname)   
    if nucleus.size == 0:
      a0 += (max_radius-a0)/Nmax; b0 = a0;
      continue

    # Find energy for frank partial with nucleus removed
    energy = setget_db(nucleus, strain, cohesv, dirname, step0)
    print("Energy of step {0} = {1}. a0 = {2}".format(step, energy, a0))
    with open(dirname+"potential_new2.dat", "a") as fp:
      fp.write(str(energy)+"\n")

    a0 += (max_radius-a0)/Nmax; b0 = a0;
    step0 += 1

#Greedy search for each new atom to add on. 
#Always choose the lowest energy state. 
def status4(strain, dirname, cohesv): 
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, tmp1, tmp2, pairs)= load_status_0(dirname) 

  Nmax = 33
  pt1pt3 = np.load(dirname+"pt1pt3.npy")
  pt1 = pt1pt3[0:3]
  a0 = 0.05; b0 = a0; max_radius = 0.75
  step0 = 0
  nucleus = select_atoms_within_ecllipse(pt1, a0, b0, dirname)   
  energy = setget_db(nucleus, strain, cohesv, dirname, step0)

  greedy_path = { }

  # Generate a series of ellipse loops, relax with strain
  for step in range(Nmax):
    bdyatoms = find_bdy_atoms(nbrlist, nucleus, totIdx, cfg)
    #assert(bdyatoms.size %2 == 0), "bdyatoms.size = {0}".format(bdyatoms.size)
    print("bdyatoms = {0}".format(bdyatoms))
    writecncfg(cfg[bdyatoms,:], H, dirname+'bdyatoms-'+str(step0))

    # Find the atom that has lowest energy if added to nucleus
    MINenergy = 1e8
    MINatom_I = -1
    MINatom_J = -1
    istep0 = 0
    bdyatoms_u = np.intersect1d(pairs[:,0], bdyatoms)
    writecncfg(cfg[nucleus,:], H, dirname+'nucleus-'+str(step0))
    with open(dirname+"B.log", "a") as fp:
      fp.write("nucleus = {0}\n".format(nucleus))
      fp.write("bdyatoms_u.size = {0}\n".format(bdyatoms_u.size))

    for atom_I in bdyatoms_u:
      istep0 -= 1
      i,j = np.where(pairs==atom_I)
      assert(j==0)
      atom_J = pairs[i[0],1]
      energy = setget_db(np.append(nucleus, [atom_I,atom_J]), strain, cohesv, dirname,'NONE' )
      writecncfg(cfg[np.append(nucleus, [atom_I,atom_J]),:], H, dirname+'nucleus-'+str(step0)+'_'+str(istep0))
      if energy < MINenergy:
        MINenergy = energy
	MINatom_I = atom_I
	MINatom_J = atom_J
        with open(dirname+"B.log", "a") as fp:
          fp.write("atom_I = {0}, atom_J = {1}, energy = {2}\n".format(atom_I, atom_J, energy))

    assert(MINatom_I >=0 and MINatom_J >=0)
    nucleus = np.append(nucleus, [MINatom_I,MINatom_J])

    with open(dirname+"B.log", "a") as fp:
      fp.write("MINatom_I = {0}, MINatom_J = {1}\n".format(MINatom_I, MINatom_J))
    with open(dirname+"potential.dat", "a") as fp:
      fp.write(str(MINenergy)+"\n")

    greedy_path[bits2str(nucleus2bits(nucleus,totIdx))] = MINenergy
  
    step0 += 1
  # end of for step in range(Nmax)
  return greedy_path

#Same as status 4 but for N atoms to add on at the same time. 
#Always choose the lowest N energy state. 
def status5(strain, dirname, cohesv): 
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, tmp1, tmp2, pairs)= load_status_0(dirname) 

  N0 = 6
  Nmax = 10000
  pt1pt3 = np.load(dirname+"pt1pt3.npy")
  pt1 = pt1pt3[0:3]
  a0 = 0.05; b0 = a0; max_radius = 0.75
  step0 = 0
  nucleus = select_atoms_within_ecllipse(pt1, a0, b0, dirname)   
  energy = setget_db(nucleus, strain, cohesv, dirname, step0)

  # Generate a series of ellipse loops, relax with strain
  for step in range(Nmax):
    bdyatoms = find_bdy_atoms(nbrlist, nucleus, totIdx, cfg)
    #assert(bdyatoms.size %2 == 0), "bdyatoms.size = {0}".format(bdyatoms.size)
    with open(dirname+"B.log", "a") as fp:
      fp.write("bdyatoms = {0}".format(bdyatoms))
    writecncfg(cfg[bdyatoms,:], H, dirname+'bdyatoms-'+str(step0))

    # Find the atoms that has lowest N energy if added to nucleus
    istep0 = 0
    bdyatoms_u = np.intersect1d(pairs[:,0], bdyatoms)
    writecncfg(cfg[nucleus,:], H, dirname+'nucleus-'+str(step0))
    MINs = -1*np.ones((bdyatoms_u.size, 3)) #energy, atom_I, atom_J
    N = min(N0, bdyatoms_u.size)
    with open(dirname+"B.log", "a") as fp:
      fp.write("nucleus = {0}\n".format(nucleus))
      fp.write("bdyatoms_u.size = {0}\n".format(bdyatoms_u.size))

    for atom_I in bdyatoms_u:
      istep0 -= 1
      i,j = np.where(pairs==atom_I)
      assert(j==0)
      atom_J = pairs[i[0],1]
      energy = setget_db(np.append(nucleus, [atom_I,atom_J]), strain, cohesv, dirname,'NONE' )
      writecncfg(cfg[np.append(nucleus, [atom_I,atom_J]),:], H, dirname+'nucleus-'+str(step0)+'_'+str(istep0))
      MINs[-istep0-1, 0] = energy
      MINs[-istep0-1, 1] = atom_I
      MINs[-istep0-1, 2] = atom_J

      with open(dirname+"B.log", "a") as fp:
        fp.write("atom_I = {0}, atom_J = {1}, energy = {2}\n".format(atom_I, atom_J, energy))

#    print MINs
    MINs = MINs[ np.argsort(MINs[:,0]), :]
#    print MINs
    nucleus = np.append(nucleus, [ MINs[0:N,1].astype(int),MINs[0:N,2].astype(int) ])

    with open(dirname+"B.log", "a") as fp:
      fp.write("MINatom_I = {0}, MINatom_J = {1}\n".format(MINs[0:N,1], MINs[0:N,2]))
    with open(dirname+"potential.dat", "a") as fp:
      fp.write(str(MINs[0:N,0])+"\n")
  
    step0 += 1

def cost(path):
  return np.max(np.array(list(path.values())))

def acceptance_probability(old_cost, new_cost, T):
  if new_cost < old_cost:
    return 1
  if new_cost >= old_cost:
    return 0.5

def neighbor(path, pairs, nbrlist, totIdx, cfg, strain, cohesv, dirname):
  new_path = path
  (N1,Y1,X1) = order_path(path)
  #choose node to perturb and form a new path
  #we don`t want to perturb the initial and end nodes
  node = 0;
  while node ==0 or node == len(N1)-1:
    node = int(math.floor(random.random()*len(N1)))
  
  nucleus = bits2nucleus(str2bits(X1[node]), totIdx)
  bdyatoms = find_bdy_atoms(nbrlist, nucleus, totIdx, cfg)
  
  bdyatoms_u = np.intersect1d(pairs[:,0], bdyatoms)
  atom_I = bdyatoms_u[int(math.floor(random.random()*len(bdyatoms_u)))]
  i,j = np.where(pairs==atom_I)
  atom_J = pairs[i[0],1]
  
  energy = setget_db(np.append(nucleus,[atom_I,atom_J]),strain,cohesv,dirname,'NONE' )
  new_path[X1[node]] = energy

  return new_path


# First call status4 and find greedy searched path as a starting point
# Then find the lowest bounds for all configurations of the same atom number
# Use simulated annealing to do global search
def status6(strain, dirname, cohesv): 
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, tmp1, tmp2, pairs)= load_status_0(dirname) 

  # find the greedy path and build data base 
  gp_file = dirname+"greedy_path";
  if os.path.exists(gp_file+".pkl"):
    greedy_path= load_obj(gp_file)
  else:
    greedy_path= status4(strain, dirname, cohesv)
    save_obj(greedy_path, gp_file)

  write_path(greedy_path, dirname+"greedy_path.log")

  # find all configurations with the same atom number 
  db_file = dirname+"db_"+strain+".2rm";
  db = load_obj(db_file)

  lowbds_record = { }
  path_record = { }
  lowbds_path = { }
  for key, val in db.iteritems():
    n = nbits1(str2bits(key))
    if n not in lowbds_record or val < lowbds_record[n]:
      lowbds_record[n] = val
      path_record[n] = key

  for n, val in lowbds_record.iteritems():
    lowbds_path[  path_record[n] ] = val

  write_path(lowbds_path, dirname+"lowbds_path.log")
  quit() 

  # Use simulated annhealing to find a global minimum
  # Starting with greedy_path as an initial condition 
  path = greedy_path
  old_cost = cost(path)
  T = 1.0
  T_min = 0.00001
  alpha = 0.9
  while T > T_min:
    i = 1
    while i <= 100:
      new_path = neighbor(path, pairs, nbrlist, totIdx, cfg, strain, cohesv, dirname)
      new_cost = cost(new_path)
      ap = acceptance_probability(old_cost, new_cost, T)
      if ap > random.random():
        path = new_path
        old_cost = new_cost
      i += 1
    T = T*alpha

  return path, cost

 
def unitest():
  status = 4
  strain = "0.05"
  dirname = '../../../runs/frankMEP/greedy_search_'+ strain +'/'
  # Change cohesive energy for different strain level

  if strain == "0.04":
    cohesv = 4.529;
  elif strain == "0.05":
    cohesv = 4.521;
  elif strain == "0.06":
    cohesv = 4.511;
  (totIdx, atomIdx, cfg, H, normal, d, pt0, nbrlist, slice_nbrlist_u, slice_nbrlist_d, pairs)= load_status_0(dirname) 
  Nmax = 200
  pt1pt3 = np.load(dirname+"pt1pt3.npy")
  pt1 = pt1pt3[0:3]
  a0 = 0.02; b0 = a0; max_radius = 0.75
  step0 = 0
  nucleus = select_atoms_within_ecllipse(pt1, a0, b0, dirname)   
  return (nucleus, nbrlist)
  
#Umbrella Sampling with an initial path
def status7():
  print("Unimplemented yet: umbrella sampling")

def main():
  strain = "0.05"
  dirname = 'runs/frankMEP/greedy_search_'+ strain +'/'
#  dirname = 'runs/frankMEP/test_db_'+ strain +'_5/' # remove two atoms
  # Change cohesive energy for different strain level

  if strain == "0.04":
    cohesv = 4.529;
  elif strain == "0.05":
    cohesv = 4.521;
  elif strain == "0.06":
    cohesv = 4.511;

  status6(strain, dirname, cohesv)
  
if __name__== "__main__":
  main()
