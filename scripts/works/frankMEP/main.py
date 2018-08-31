import os, sys
import numpy as np
import pickle
import math
import random
from subprocess import call
import numpy as np
from numpy import linalg as LA
from utility import *

from MDobj import MDobj
from SearchAlgorithms import GreedySearch
from SearchAlgorithms import RlSearch

sys.path.append('../../lib/')
sys.path.append('../../tests/')
import mdsw
from View import Viewer

# Control parameters we need to set for Frank Partial Dislocation MEP calculation

strain = "0.05"
dirname = 'greedy_search_'+ strain +'/'

# Cohesive energy is determined by applied strain.
# Run separate simualations to calcualte these parameters.

if strain == "0.04":
    cohesv = 4.529;
elif strain == "0.05":
    cohesv = 4.521;
elif strain == "0.06":
    cohesv = 4.511;

# Maximum iterations for Greedy Search algorithms.

Nmax = 33

# Color codes for View module

red =   [1.0, 0.0, 0.0, 1.0]
green = [0.0, 1.0, 0.0, 1.0]
blue =  [0.0, 0.0, 1.0, 1.0]

# Initialize alloc all variables

swobj = MDobj(strain, dirname, cohesv)
swobj.initialize()
swobj.sw.writeatomeyecfg("tf.cfg")

# Choose which search algorithm to use

alg = GreedySearch(cohesv, strain, dirname)
#alg = RlSaerch(swobj, cohesv, strain, dirname)

# A search begins.

alg.search(swobj, Nmax)

#alg.visualize_path()


if 0:
    nucleus = swobj.choose_elipse_state(np.array([0, 0, 0.38475]), 0.35, 0.35)
    plotlist =np.extract(np.abs(swobj.SR[:,2])>0.375, np.arange(swobj.sw.NP))

if 0:
    pltlist = np.concatenate((plotlist,swobj.slice_nbrlist_d, swobj.slice_nbrlist_u))
    colorlist = np.vstack((np.tile(red,(len(plotlist),1)), np.tile(blue,(len(swobj.slice_nbrlist_d),1)), np.tile(green,(len(swobj.slice_nbrlist_u),1))))

if 0:
    pltlist = np.concatenate((plotlist,swobj.slice_nbrlist_d, nucleus))
    colorlist = np.vstack((np.tile(red,(len(plotlist),1)), np.tile(blue,(len(swobj.slice_nbrlist_d),1)), np.tile(green,(len(nucleus),1))))

if 0:
    pltlist = np.concatenate((plotlist,swobj.slice_nbrlist_u, nucleus))
    colorlist = np.vstack((np.tile(red,(len(plotlist),1)), np.tile(blue,(len(swobj.slice_nbrlist_u),1)), np.tile(green,(len(nucleus),1))))


if 0:
    pltlist = np.concatenate((plotlist,idx_u, nucleus))
    colorlist = np.vstack((np.tile(red,(len(plotlist),1)), np.tile(blue,(len(idx_u),1)), np.tile(green,(len(nucleus),1))))

if 0:
    pltlist = np.concatenate((plotlist,idx_d, nucleus))
    colorlist = np.vstack((np.tile(red,(len(plotlist),1)), np.tile(blue,(len(idx_d),1)), np.tile(green,(len(nucleus),1))))

if 0:
    pltlist = np.concatenate((plotlist,idx_u, idx_d))
    colorlist = np.vstack((np.tile(red,(len(plotlist),1)), np.tile(blue,(len(idx_u),1)), np.tile(green,(len(idx_d),1))))
    view = Viewer(swobj, 300, 300, pltlist, colorlist)
    view.rendering()

if 0:
    red =   [1.0, 0.0, 0.0, 1.0]
    green = [0.0, 1.0, 0.0, 1.0]
    blue =  [0.0, 0.0, 1.0, 1.0]
    surfnds=np.extract(np.abs(self.SR[:,2])>0.375, np.arange(self.sw.NP))
    #fixednodes=fixednodes_d #np.extract(self.fixed==1, np.arange(self.sw.NP))
    group = self.sw.group()
    group1nodes = np.extract(group==1, np.arange(self.sw.NP))
    plotlist = np.concatenate((surfnds, group1nodes))
    colorlist = np.vstack((np.tile(red,(len(surfnds),1)), np.tile(blue,(len(group1nodes),1))))
    view = Viewer(self, 300, 300, plotlist, colorlist)
    view.rendering()
    self.sw.sleep()



#swobj.step(nucleus)
#swobj.step(nucleus)
#swobj.step(nucleus)



#view = Viewer(swobj, 300, 300, np.union1d(np.union1d(plotlist,swobj.slice_nbrlist_d), swobj.slice_nbrlist_u))
#view = Viewer(swobj, 300, 300, np.union1d(np.union1d(np.union1d(plotlist,nucleus), idx_u), idx_d))
#view = Viewer(swobj, 300, 300) , np.union1d(np.union1d(plotlist,idx_u), idx_d))
#exit(0)


if 0:
    print(swobj.totIdx.shape)
    print(nucleus.shape)

    swobj1 = MDobj(strain, dirname, cohesv)
    swobj1.reset()
    print("1")
    print(swobj1.SR.shape)
    print("2")
    tmp = np.copy(swobj1.SR)
    print(tmp.shape)
    #swobj1.SR = np.copy(swobj.SR[swobj.totIdx,:])
    #np.copyto(swobj1.SR, swobj.SR[swobj.totIdx,:])
    np.copyto(swobj1.SR, tmp)
    print("3")
    print(swobj1.SR.shape)
    print("4")
    swobj1.sw.writeatomeyecfg("swobj1.cfg")

    plotlist =np.extract(np.abs(swobj.SR[:,2])>0.35, np.arange(swobj.SR.shape[0]))
    #view = Viewer(swobj, 600, 600, swobj.totIdx )
    view = Viewer(swobj, 600, 600, np.union1d(plotlist, nucleus) )
    view.rendering()
    swobj1.sw.sleep()

print("Program exits normally")


