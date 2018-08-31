""" MDobj is a wrapper of mdsw + a set of nbrlist operations.
"""
import os
import sys
sys.path.append('../../../lib/')
sys.path.append('../../tests/')
import mdsw
from numpy import linalg as LA
from View import Viewer
from utility import *

class MDobj(object):
    def __init__(self, strain, dirname, cohesv):
        self.strain = np.double(strain)
        self.dirname = dirname
        self.sw = mdsw.SWFrame()
        self.SR = np.array([])
        self.fixed = np.array([])
        self.group = np.array([])
        self.slice_nbrlist_u = np.array([])
        self.slice_nbrlist_d = np.array([])
        self.totIdx = np.array([])
        self.nbrlist = np.array([])
        self.nn = np.array([])
        self.normal = np.array([0, 0, 0])
        self.cohesv = cohesv
        self.pairs = np.array([])
        self.E0 = 0 # store E of perfect structure

    def relax_fixbox(self):
        self.sw.conj_ftol = 1e-4
        self.sw.conj_itmax = 3800
        self.sw.conj_fevalmax =2 # 6000
        self.sw.conj_fixbox = 1
        self.sw.relax()

    """ Find atoms eps away from a plane defined by its normal (nrml) and a _pt
        If opt > 0, choose atoms above the plane
        If opt < 0, choose atoms below the plane
        If opt = 0, returns them all
    """
    def select_totalatoms_on_slice(self, eps, _pt, nrml, opt):
        idx = np.extract((np.abs(nrml[0]*self.SR[:,0]+nrml[1]*self.SR[:,1]+nrml[2]*self.SR[:,2]-np.ones(self.sw.NP)*np.dot(_pt,nrml))<eps),np.arange(self.sw.NP))
        if opt == 0:
            return idx
        else:
            D = np.ones(len(idx))*np.dot(_pt,nrml)
        cond=nrml[0]*self.SR[idx,0]+nrml[1]*self.SR[idx,1]+nrml[2]*self.SR[idx,2]-D>0 if opt>0 else nrml[0]*self.SR[idx,0]+nrml[1]*self.SR[idx,1]+nrml[2]*self.SR[idx,2]-D<0
        return np.extract(cond,idx)


    # Prepare surface reconstructed thin film geometry
    def initialize(self):

        # Create a perfect si lattice
        self.sw.initvars()
        self.sw.dirname = self.dirname
        self.sw.NNM = 100
        self.sw.crystalstructure = "diamond-cubic"
        self.sw.latticeconst = mdsw.VectorDouble([5.4309529817532409, 5.4309529817532409, 5.4309529817532409])
        # <e_x^1,e_x^2,e_x^3,nx>
        self.sw.latticesize = mdsw.VectorDouble([ 1 ,1, 0, 12,  -1, 1, 0, 12,  0, 0, 1, 16 ])
        self.sw.makecrystal()
        self.sw.saveH()
        self.relax_fixbox()

        # Create a free surface by making vaccumm region
        self.sw.vacuumratio = 1.0-1.0/(1.0+0.2)
        self.sw.conj_fixbox = 0
        self.sw.input = mdsw.VectorDouble([ 3, 3, 0.2 ])
        self.sw.changeH_keepR()
        self.relax_fixbox()

        print("initialize I am here 0")
        # Surface reconstruction in [110] directions
        ny = self.sw.latticesize[7]
        for i in range(int(ny)):
             ymin = -0.5006+1.0/ny*i
             ymax = -0.5006+1.0/ny*(i+0.5)
             self.sw.input = mdsw.VectorDouble([ 1,-10,10,ymin,ymax,0.403,10 ])
             self.sw.fixatoms_by_position()
             self.sw.input = mdsw.VectorDouble([ 1,-10,10,ymin,ymax,-10,-0.416 ])
             self.sw.fixatoms_by_position()

        self.sw.input = mdsw.VectorDouble([1])
        self.sw.setfixedatomsgroup()
        self.sw.freeallatoms()

        for i in range(int(ny)):
             ymin = -0.5006+1.0/ny*(i+0.5)
             ymax = -0.5006+1.0/ny*(i+1)
             self.sw.input = mdsw.VectorDouble([ 1,-10,10,ymin,ymax,0.403,10 ])
             self.sw.fixatoms_by_position()
             self.sw.input = mdsw.VectorDouble([ 1,-10,10,ymin,ymax,-10,-0.416 ])
             self.sw.fixatoms_by_position()

        self.sw.input = mdsw.VectorDouble([2])
        self.sw.setfixedatomsgroup()
        self.sw.freeallatoms()
        print("initialize I am here 1")

        self.sw.input = mdsw.VectorDouble([ 1,  0,  0.8, 0,  1 ])
        self.sw.movegroup()
        self.sw.input = mdsw.VectorDouble([ 1,  0, -0.8, 0,  2 ])
        self.sw.movegroup()
        print("initialize I am here 2")
        self.group.fill(0)
        print("initialize I am here 3")

        # Relax and save perfect thin film configurations
        self.relax_fixbox()

        # Make a reference of commonly used arrays in sw
        self.SR = self.sw.SR()
        self.fixed = self.sw.fixed()
        self.group = self.sw.group()
        print("initialize I am here 4")

        # Make copies of initial states to go back to anytime
        self.sw.setconfig1()
        self.sw.saveH()
        self.NP0 = self.sw.NP
        self.sw.eval()
        self.E0 = self.sw.EPOT

        # Prepare nbrlist of slip planes
        self.make_nnlist2d()

    """ Prepare nn list for neighbor atom ids of each atom on the 2d plane
        pick the contrl pts to define a slice of atoms
        atoms between left and right glide plane will be removed.
        pt1_u/pt1_d----pt2_u/pt2_d  --> left glide plane
        pt3_u/pt3_d----pt4_u/pt4_d  --> right glide plane
        pt0: middle of pt1 and pt3, thus on shuffle plane
    """
    def make_nnlist2d(self):
        pt1d = np.array([0.0208,0, 0.3784]); pt1u = np.array([0.0208,0, 0.3911])
        pt2d = np.array([0.1042, 0, 0.2734]);pt2u = np.array([0.1042, 0, 0.2865])
        pt3d = np.array([0.0625,0, 0.3773]); pt3u = np.array([0.0625,0, 0.3911])
        pt4d = np.array([0.1250, 0, 0.2995]);pt4u = np.array([0.1250, 0, 0.3125])
        pt1 = 0.5*(pt1d + pt1u);  pt2 = 0.5*(pt2d + pt2u)
        pt3 = 0.5*(pt3d + pt3u);  pt4 = 0.5*(pt4d + pt4u)
        pt0 = 0.5*(pt1+pt3);
        print("pt1 = {0}, pt3 = {1}".format(pt1, pt3))

        v = (pt2-pt1).astype(np.double)
        v /=LA.norm(v)
        self.normal = np.cross(v, [0,1,0])

        eps = 0.02
        totIdx_1 = self.select_totalatoms_on_slice(eps, pt1, self.normal, 1)
        totIdx_2 = self.select_totalatoms_on_slice(eps, pt3, self.normal, -1)
        self.totIdx = np.union1d(totIdx_1, totIdx_2)
        self.pairs  = self.make_pairs(totIdx_1, totIdx_2)

        # Build neighborhood list for totIdx
        self.sw.refreshnnlist()
        # self.sw.fprintnnlist()
        self.nbrlist = self.sw.nindex() # this actually bind to nindex_mem
        self.nn = self.sw.nn()

        # Select atoms sandwich pt0-normal plane with 2*eps thickness
        totIdx_u = self.select_totalatoms_on_slice(2*eps, pt0, self.normal, 1)
        totIdx_d = self.select_totalatoms_on_slice(2*eps, pt0, self.normal, -1)

        # Select atoms above/below pt0-normal plane minus nucleus atoms
        # Select atoms to be moved when closing-up a trench is subset of these
        self.slice_nbrlist_u = np.intersect1d(getnbrlist(totIdx_1, self.nbrlist), totIdx_u)
        self.slice_nbrlist_d = np.intersect1d(getnbrlist(totIdx_2, self.nbrlist), totIdx_d)

    """ Find closest pair of atoms across the glide-set plane
    """
    def make_pairs(self, totIdx_1, totIdx_2):
        (t1,t2) = (totIdx_1,totIdx_2) if len(totIdx_1)<len(totIdx_2) else (totIdx_2,totIdx_1)
        pairs_ = [ ]
        for atom_I in t1:
            mindis = 1e5
            for atom_J in t2:
                tmp = LA.norm(self.SR[atom_I,:]-self.SR[atom_J,:])
                if tmp < mindis:
                    mindis = tmp
                    minatom = atom_J
            pairs_.append([atom_I, minatom])
        return np.array(pairs_).astype(int)

    def make_frk_dislocation(self, nucleus):
        idx_u = np.intersect1d(getnbrlist(nucleus, self.nbrlist), self.slice_nbrlist_u)
        idx_d = np.intersect1d(getnbrlist(nucleus, self.nbrlist), self.slice_nbrlist_d)

        # Perturb atoms on both sides of nucleus (within nbrlist)
        print("frk I am here 1")
        self.sw.freeallatoms()
        for id in [ idx_u ]: self.fixed[id] = 1
        self.sw.input = mdsw.VectorDouble([1])
        self.sw.setfixedatomsgroup()
        self.sw.freeallatoms()
        print("frk I am here 2")
        for id in [ idx_d ]: self.fixed[id] = 1
        self.sw.input = mdsw.VectorDouble([2])
        self.sw.setfixedatomsgroup()
        self.sw.freeallatoms()
        print("frk I am here 3")

        mag  = 0.9
        magx = self.normal[0] * mag
        magy = self.normal[1] * mag
        magz = self.normal[2] * mag

        print("frk I am here 4")
        self.sw.input = mdsw.VectorDouble([ 1, -magx, -magy, -magz, 1 ])
        self.sw.movegroup()
        self.sw.input = mdsw.VectorDouble([ 1, magx, magy, magz, 2 ])
        self.sw.movegroup()

        print("frk I am here 5")
        # Put group back to 0. Important!
        self.group.fill(0)

        # Remove nucleus from SR
        for id in nucleus: self.fixed[id] = 1
        self.sw.removefixedatoms()
        self.sw.freeallatoms()

        print("frk I am here 6")
        # Apply strain to close trench and create a frank partial
        H = self.sw.H
        H[0] = self.sw.H[0]*(1.0-self.strain)
        self.sw.H = H
        self.relax_fixbox()
        self.sw.SHtoR()
        print("frk I am here 7")

        # Evaluate and return potential energy
        self.sw.eval()
        return self.sw.EPOT

    """ Select atoms belong to an eclipse defined by a, b, pt0 and normal
    """
    def choose_elipse_state(self, pt0, a, b):
        e0 = [1,0,0]; e1 = [0,1,0]; e2 = [0,0,1];
        ep0 = e1;     ep1 = np.cross(self.normal, ep0); ep2 = self.normal;

        # Coordinate transform
        Q = np.zeros((3,3))
        Q[0][0]=np.dot(ep0,e0);Q[1][1]=np.dot(ep1,e1);Q[2][2]=np.dot(ep2,e2)
        Q[0][1]=np.dot(ep0,e1);Q[0][2]=np.dot(ep0,e2);Q[1][0]=np.dot(ep1,e0)
        Q[1][2]=np.dot(ep1,e2);Q[2][0]=np.dot(ep2,e0);Q[2][1]=np.dot(ep2,e1)

        return np.intersect1d(np.extract(np.dot(Q, (self.SR-pt0).transpose())[0,:]**2/a**2+np.dot(Q, (self.SR-pt0).transpose())[1,:]**2/b**2<=1, np.arange(self.sw.NP)), self.totIdx)

    # Reset to have the initial state
    def reset(self):
        # Make these random later: choose an initial small nucleus
        pt0 = np.array([0,0,0])
        a0 = 0.05
        b0 = a0
        nucleus = self.choose_elipse_state(pt0, a0, b0)
        return nucleus


    """ step() returns four values
        observation: state of environemnt
        reward of previous action.
        done: True if the episode is finished
        info: diagnostic information for debug
        Retrieve energy from a database if the state (bitarray) is already there
        Otherwise, calculate the state and store it to the data base
    """

    def step(self, nucleus):

        db_file = self.dirname+"db_" + str(self.strain);
        db  = load_obj(db_file) if os.path.exists(db_file+".pkl") else { }
        bitstr = bits2str(nucleus2bits(nucleus, self.totIdx))

        if bitstr in db:
            print("data base has = {0} data points".format(len(db)))
            return nucleus, db[bitstr], db[bitstr]<self.E0, {}
        else:
            # Perturb atoms on boths sides of nucleus
            print("step I am here 0")

            energy = self.make_frk_dislocation(nucleus)
            energy -= self.cohesv * nucleus.size;
            print("step I am here 1")
            db[bitstr] = energy
            save_obj(db, db_file)
            print("step I am here 2")

            # Put SR, H back
            self.sw.NP = self.NP0
            self.sw.SR1toSR()
            print("step I am here 3")
            self.sw.restoreH()
            print("step I am here 4")
            self.sw.refreshnnlist()
            print("step I am here 5")

        # If energy drops below initial energy, done!
        return nucleus, energy, energy<self.E0, {}


