""" Search Algorihtms: 1) Greedy search. 2) Deep Q-learning.
"""
from MDobj import MDobj
from utility import *

class GreedySearch(object):
    def __init__(self, swobj, cohesv, strain, dirname): 
        self.cohesv = cohesv
        self.strain = strain
        self.dirname = dirname
        self.swobj = swobj

    #Greedy search for each new atom to add on: always choose the lowest energy state. 
    def search(self, Nmax): 
        pt1 = self.pt1pt3[0:3]
        step0 = 0
        nucleus = swobj.reset()
        s, energy, _,_ = self.step(nucleus)
    
        self.MEP = { }
    
        # Generate a series of ellipse loops, relax with strain
        for step in range(Nmax):
            bdyatoms = find_nbr_atoms(self.swobj.nbrlist, self.sw.state, self.sw.totIdx)
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
                nucleus = np.array(np.append(nucleus, [atom_I, atom_J]))
                nucleus, energy, done, info = self.swobj.step(nucleus)
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
    
            self.MEP[bits2str(nucleus2bits(nucleus,totIdx))] = MINenergy
        
            step0 += 1 # end of for step in range(Nmax)
    
        return self.MEP
    

    #def visualize_path():
    #    (X,Y,Z) = order_path(self.MEP)
    #    plt.plot(X,Y,'*-')
    #    plt.draw()
    
    def write_path(path, name):
        (N,Y,X) = order_path(self.MEP)
        with open(name, "w") as fp:
            for val in N:
                fp.write(str(val)+" ")
            fp.write("\n")
            for val in Y:
                fp.write(str(val)+" ")
            fp.write("\n")
      
    def order_path():
        X = np.array(list(self.MEP.keys()))
        Y = np.array(list(self.MEP.values()))
        N = [ ]
        for item in X:
            N.append(nbits1(str2bits(item)))
        N = np.array(N) 
        inds = N.argsort()
        Y = Y[inds]
        N = N[inds]
        X = X[inds]
        return (N, Y, X) 


class RlSearch(object):
    def __init__(self, swobj, cohesv, strain, dirname): 
        self.cohesv = cohesv
        self.strain = strain
        self.dirname = dirname
        self.swobj = swobj

    def search(self, Nmax):
        raise NotImplementedError

