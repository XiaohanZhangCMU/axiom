import os, sys
import numpy as np
import pickle
import math
import random
from subprocess import call
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

# Automatically reuse atom container if available. Rewrite if this is not the intention.

swobj = MDobj(strain, dirname, cohesv)
swobj.reset()
nucleus = swobj.choose_elipse_state(np.array([0, 0, 0.38475]), 0.32, 0.32)

swobj.step(nucleus)
swobj.step(nucleus)
swobj.step(nucleus)

swobj.sw.writeatomeyecfg("test_step_done.cfg")
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

# A set of test functions for MDobj class

alg = GreedySearch(cohesv, strain, dirname)
alg.search(swobj, Nmax)

#alg = RlSaerch(swobj, cohesv, strain, dirname)

# A set of test functions for Reinforcement Learning algorithm


# A search begins.

#alg.search()

#alg.visualize_path()
print("Program exits normally")


