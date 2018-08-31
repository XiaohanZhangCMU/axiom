""" Search Algorihtms: 1) Greedy search. 2) Deep Q-learning.
"""
from MDobj import MDobj
from utility import *
from View import Viewer
from MDobj import MDobj

class GreedySearch(object):
    def __init__(self, cohesv, strain, dirname):
        self.cohesv = cohesv
        self.strain = strain
        self.dirname = dirname

    #Greedy search for each new atom to add on: always choose the lowest energy state.
    def search(self, swobj, Nmax):
        step0 = 0
        saveinter = 0

        nucleus = swobj.choose_elipse_state(np.array([0, 0, 0.38475]), 0.08, 0.08)
        s, energy, _,_ = swobj.step(nucleus)

        self.MEP = { }

        # Generate a series of ellipse loops, relax with strain
        for step in range(Nmax):
            bdyatoms = find_nbr_atoms(swobj.nbrlist, nucleus, swobj.totIdx)
            #assert(bdyatoms.size %2 == 0), "bdyatoms.size = {0}".format(bdyatoms.size)
            #print("bdyatoms = {0}".format(bdyatoms))

            # Find the atom that has lowest energy if added to nucleus
            MINenergy = 1e8
            MINatom_I = -1
            MINatom_J = -1
            istep0 = 0
            bdyatoms_u = np.intersect1d(swobj.pairs[:,0], bdyatoms)

            nucleus0 = np.copy(nucleus)

            for atom_I in bdyatoms_u:
                istep0 -= 1
                i,j = np.where(swobj.pairs==atom_I)
                assert(j==0)
                atom_J = swobj.pairs[i[0],1]
                nucleus = np.append(nucleus, [atom_I, atom_J])
                nucleus, energy, done, info = swobj.step(nucleus)

                if energy < MINenergy:
                    MINenergy = energy

                MINatom_I = atom_I
                MINatom_J = atom_J

                with open(swobj.dirname+"B.log", "a") as fp:
                    fp.write("nucleus size = {0}; (atom I, atom J) = ({1}); energy = {2}\n".format(len(nucleus), (atom_I, atom_J), energy))

            with open(swobj.dirname+"B.log", "a") as fp:
                fp.write("------------------------------\n")

            assert(MINatom_I >=0 and MINatom_J >=0)
            nucleus = np.append(nucleus0, [MINatom_I,MINatom_J])

            with open(swobj.dirname+"potential.dat", "a") as fp:
                fp.write(str(MINenergy)+"\n")
            if 1:
                self.save_path_node(swobj, nucleus, saveinter)
                saveinter+=1

            self.MEP[bits2str(nucleus2bits(nucleus,swobj.totIdx))] = MINenergy

            save_obj(self.MEP, 'mep')

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

    def visualize_path_node(self, swobj, nucleus):
        red =   [1.0, 0.0, 0.0, 1.0]
        green = [0.0, 1.0, 0.0, 1.0]
        blue =  [0.0, 0.0, 1.0, 1.0]
        plotlist =np.extract(np.abs(swobj.SR[:,2])>0.375, np.arange(swobj.sw.NP))
        atomlist = np.concatenate((plotlist,nucleus))
        colorlist = np.vstack((np.tile(red,(len(plotlist),1)), np.tile(blue,(len(nucleus),1))))
        view = Viewer(swobj, 300, 300, atomlist, colorlist)
        view.rendering()

    def save_path_node(self, swobj, nucleus, saveinter):
        swobj.sw.finalcnfile=swobj.dirname + "/pathimg_"+str(saveinter)+".cfg"
        swobj.sw.freeallatoms()
        swobj.fixed.fill(1)
        swobj.fixed[nucleus] = 0
        swobj.sw.removefixedatoms()
        swobj.sw.writeatomeyecfg(swobj.sw.finalcnfile)

        swobj.sw.NP = swobj.NP0
        swobj.sw.SR1toSR()
        swobj.sw.refreshnnlist()


class DQNSearch(object):
    def __init__(self, swobj, cohesv, strain, dirname):
        self.cohesv = cohesv
        self.strain = strain
        self.dirname = dirname

    def search(self, Nmax):
        raise NotImplementedError

