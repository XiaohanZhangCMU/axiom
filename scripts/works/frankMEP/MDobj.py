""" A class to hold geometry of atoms through mdsw + a set of functions to select atoms and build neighborhoods.
"""
import sys
sys.path.append('../../../lib/')
sys.path.append('../../tests/')
import mdsw
from View import Viewer
from utility import *

class MDobj(object):
    def __init__(self, strain, dirname):
        self.strain = strain
        self.dirname = dirname
        self.sw = mdsw.SWFrame()
        self.SR0 = np.array([])
        self.slice_nbrlist_u = np.array([])
        self.slice_nbrlist_d = np.array([])
        self.nbrlist = np.array([])
        self.nn = np.array([])

    
    def relax_fixbox(self):
        self.sw.conj_ftol = 1e-4
        self.sw.conj_itmax = 3800
        self.sw.conj_fevalmax = 6000
        self.sw.conj_fixbox = 1
        self.sw.relax()

    # Find atoms eps away from a plane defined by its normal (nrml) and a _pt
    # If opt > 0, choose atoms above the plane
    # If opt < 0, choose atoms below the plane
    # If opt = 0, returns them all

    def select_totalatoms_on_slice(self, eps, _pt, normal, opt):
        idx = np.extract((np.abs(nrml[0]*self.SR[:,0]+nrml[1]*self.SR[:,1]+nrml[2]*self.SR[:,2]-np.ones(self.SR.shape[0])*np.dot(_pt,nrml))<eps),np.xarrange(self.NP0))
        if opt == 0: 
            return idx 
        else:
            D = np.ones(len(idx))*np.dot(_pt,nrml)
        cond=nrml[0]*self.SR[idx,0]+nrml[1]*self.SR[idx,1]+nrml[2]*self.SR[idx,2]-D>0 if opt>0 else nrml[0]*self.SR[idx,0]+nrml[1]*self.SR[idx,1]+nrml[2]*self.SR[idx,2]-D<0
        return np.extract(cond,idx)

    def choose_elipse_state(self, pt0, a, b):
        print("Select nucleus within an elipse.")
        e0 = [1,0,0]; e1 = [0,1,0];                e2 = [0,0,1];
        ep0 = e1;     ep1 = np.cross(normal, ep0); ep2 = normal; 

        # coordinate transform
        Q = np.zeros((3,3))
        Q[0][0] = np.dot(ep0,e0); Q[1][1] = np.dot(ep1,e1); Q[2][2] = np.dot(ep2,e2)
        Q[0][1] = np.dot(ep0,e1); Q[0][2] = np.dot(ep0,e2); Q[1][0] = np.dot(ep1,e0)
        Q[1][2] = np.dot(ep1,e2); Q[2][0] = np.dot(ep2,e0); Q[2][1] = np.dot(ep2,e1)
    
        self.state = []; 
        for i in range(self.NP0):
            _x = np.dot(Q,SR[i,:]-pt0);# origin shifts to the picked point and project 
            x = _x[0]; y = _x[1]; z = _x[2]; 
            if (x)*(x)/a/a + (y)*(y)/b/b <= 1 :
                self.sate.append(i);
        self.state = np.array(self.state)
        return np.intersect1d(self.state, totIdx)
  

    # Reset to have the initial state
    def reset(self):
        if self.SR0.size > 0: # Reuse SR0
            self.SR = np.copy(self.SR0)

        # Make these random later: choose an initial small nucleus 
        pt0 = np.array([0,0,0])
        a0 = 0.05
        b0 = a0

        # Create a perfect si lattice
        self.sw.initvars()
        self.sw.dirname = self.dirname
        self.sw.NNM = 200
        self.sw.crystalstructure = "diamond-cubic"
        self.sw.latticeconst = mdsw.VectorDouble([5.4309529817532409, 5.4309529817532409, 5.4309529817532409])
        # <e_x^1,e_x^2,e_x^3,nx>
        self.sw.latticesize = mdsw.VectorDouble([ 1 ,1, 0, 12,  -1, 1, 0, 12,  0, 0, 1, 16 ])
        self.sw.makecrystal()
        self.sw.saveH()
        self.sw.relax()

        # Create a free surface by making vaccumm region
        self.sw.vacuumratio = 1.0-1.0/(1.0+0.2)
        self.sw.conj_fixbox = 0
        self.sw.input = mdsw.VectorDouble([ 3, 3, 0.2 ])
        self.sw.changeH_keepR()
        self.sw.relax()

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

        self.sw.input = mdsw.VectorDouble([ 1,  0,  0.8, 0,  1 ])
        self.sw.movegroup()
        self.sw.input = mdsw.VectorDouble([ 1,  0, -0.8, 0,  2 ])
        self.sw.movegroup()
        
        # Relax and save perfect thin film configurations
        self.relax_fixbox()
        self.sw.finalcnfile = "0K_0.0_relaxed_surf001.cn"
        self.sw.writecn(0,False)
        self.sw.finalcnfile = "0K_0.0_relaxed_surf001.cfg"
        self.sw.writeatomeyecfg(self.sw.finalcnfile)

        # Make a deep copy of the perfect thin film configuration to SR0
        self.SR0 = np.copy(self.sw.SR())
        #self.H0  = np.copy(self.sw.H())
        self.NP0 = self.sw.NP
        
    # Prepare nn list for neighbor atom ids of each atom on the 2d plane
    # pick the contrl pts to define a slice of atoms 
    # atoms between left and right glide plane will be removed.
    # pt1_u/pt1_d----pt2_u/pt2_d  --> left glide plane
    # pt3_u/pt3_d----pt4_u/pt4_d  --> right glide plane
    # pt0: middle of pt1 and pt3, thus on shuffle plane
    def nnlist2d(self):
        pt1d = np.array([0.0208,0, 0.3784]); pt1u = np.array([0.0208,0, 0.3911])
        pt2d = np.array([0.1042, 0, 0.2734]);pt2u = np.array([0.1042, 0, 0.2865])
        pt3d = np.array([0.0625,0, 0.3773]); pt3u = np.array([0.0625,0, 0.3911])
        pt4d = np.array([0.1250, 0, 0.2995]);pt4u = np.array([0.1250, 0, 0.3125])
        pt1 = 0.5*(pt1d + pt1u);  pt2 = 0.5*(pt2d + pt2u)
        pt3 = 0.5*(pt3d + pt3u);  pt4 = 0.5*(pt4d + pt4u)
        pt0 = 0.5*(pt1+pt3); 
    
        v = (pt2-pt1)/LA.norm(v)
        normal = np.cross(v, [0,1,0])
    
        eps = 0.02
        totIdx_1 = select_totalatoms_on_slice(eps, pt1, normal, 1)
        totIdx_2 = select_totalatoms_on_slice(eps, pt3, normal, -1)
        print("pt1 = {0}, pt3 = {1}".format(pt1, pt3))
        self.totIdx = np.union1d(totIdx_1, totIdx_2)
        self.pairs  = make_pairs(totIdx_1, totIdx_2)
    
        # Build neighborhood list for totIdx
        self.nbrlist = sw.nindex() # ?????????????????/
        self.nn = sw.nn()

        # atoms sandwich pt0-normal plane with 3*eps thickness
        totIdx_u = select_totalatoms_on_slice(3*eps, pt0, normal, 1)
        totIdx_d = select_totalatoms_on_slice(3*eps, pt0, normal, -1)

        # atoms above/below pt0-normal plane minus nucleus atoms
        # atoms to be moved when closing-up a trench is subset of these
        self.slice_nbrlist_u = np.intersect1d(getnbrlist(totIdx_1, nbrlist), totIdx_u)
        self.slice_nbrlist_d = np.intersect1d(getnbrlist(totIdx_2, nbrlist), totIdx_d)
    
#        self.nbr_above = np.intersect1d(getnbrlist(state, nbrlist), slice_nbrlist_u)
#        self.nbr_below = np.intersect1d(getnbrlist(state, nbrlist), slice_nbrlist_d)
      
    
    # Then relax and compute energy which is written to EPOT_2.dat
    def make_frk_dislocation(self, nucleus):
        tol = 0.05
        SR = self.sw.SR()
        fixed = self.sw.fixed()

        # Perturb atoms on both sides of nucleus (within nbrlist) 
        for id in [ self.slice_nbrlist_u ]: fixed[id] = 1
        self.input = mdsw.VectorDouble([1])
        self.sw.setfixedatomsgroup()
        self.sw.freeallatoms()

        for id in [ self.slice_nbrlist_d ]: fixed[id] = 1
        self.sw.input = mdsw.VectorDouble([2])
        self.sw.setfixedatomsgroup()
        self.sw.freeallatoms()

        mag  = 0.9
        magx = self.normal[0] * mag 
        magy = self.normal[1] * mag 
        magz = self.normal[2] * mag 

        self.sw.input = mdsw.VectorDouble([ 1, -magx, -magy, -magz, 1 ])
        self.sw.movegroup()
        self.sw.input = mdsw.VectorDouble([ 1, magx, magy, magz, 1 ])
        self.sw.movegroup()

        # Remove nucleus from SR
        for id in self.nucleus: self.sw.fixed[id] = 1
        self.sw.removefixedatoms()
        self.sw.freeallatoms()

        # Apply strain to close trench and create a frank partial
        H11_fix = H0[0][0]*(1.0-self.strain)
        H = self.sw.H()
        H[0][0] = H11_fix #by reference, therefore sw.H has the same memory address as H
        self.relax_fixbox()
        self.sw.SHtoR()

        # Evaluate and return potential energy
        sw.eval()
        return sw.EPOT 
    
    """ step() returns four values
        observation: state of environemnt
        reward of previous action.
        done: True if the episode is finished
        info: diagnostic information for debug
        Retrieve energy from a database if the state (bitarray) is already there
        Otherwise, calculate the state and store it to the data base
    """
    
    def step(self, nucleus):

        db_file = sw.dirname+"db_" + sw.strain;
        db  = load_obj(db_file) if os.path.exists(db_file+".pkl") else { }
      
        bitstr = bits2str(nucleus2bits(nucleus, self.totIdx))
        if bitstr in db:
            print("data base has = {0} data points".format(len(db)))
            return nucleus, db[bitstr], False, {} #want to maximize the energy cost
        else:
            # Perturb atoms on boths sides of nucleus
            # Then evaluate the potential energy of the perturbed system
            energy = self.make_frk_dislocation(nucleus)
            energy -= cohesv * nucleus.size;
            db[bitstr] = energy
            save_obj(db, db_file)

            # Put back SR, H to SR0, H0 after make_frk_dislocation is called
            self.sw.SR = np.copy(SR0)
            self.sw.H = np.copy(H)

        return nucleus, energy, False, {} #want to maximize the energy cost

