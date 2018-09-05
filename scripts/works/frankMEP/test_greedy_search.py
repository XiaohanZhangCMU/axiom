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
from SearchAlgorithms import DQNSearch

sys.path.append('../../lib/')
sys.path.append('../../tests/')
import mdsw
from View import Viewer

""" Specify a small B state, check if GreedySearch works as expected.
    Run this script. Check EPOT.dat and path images.
"""

# Control parameters

strain = "0.05"
dirname = 'greedy_search_'+ strain +'/'

# Cohesive energy is determined by applied strain.
# Run separate simualations to determined Cohesive energy.

if strain == "0.04":
    cohesv = 4.529;
elif strain == "0.05":
    cohesv = 4.521;
elif strain == "0.06":
    cohesv = 4.511;

# Color codes for View module

red =   [1.0, 0.0, 0.0, 1.0]
green = [0.0, 1.0, 0.0, 1.0]
blue =  [0.0, 0.0, 1.0, 1.0]

# Initialize alloc all variables

swobj = MDobj(strain, dirname, cohesv)
swobj.initialize()

# Create State A and B

stateA = np.array([])
stateB = swobj.choose_elipse_state(np.array([0,0,0.38475]), 0.2, 0.2)

swobj.make_frk_dislocation(stateB)
energy = swobj.eval(stateB)
swobj.sw.writeatomeyecfg(dirname + '/img_B.cfg')
swobj.saveNucleus(stateB, dirname + '/nucleus_B.cfg')

swobj.restoreConfig()

# Choose which search algorithm to use

alg = GreedySearch(stateA, stateB)

# Begin search
alg.search(swobj)

print("Program exits normally")


