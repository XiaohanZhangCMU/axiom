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


""" Select an eclipse of atoms and check geometry.
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

# Maximum iterations for Greedy Search algorithms.

# Color codes for View module

red =   [1.0, 0.0, 0.0, 1.0]
green = [0.0, 1.0, 0.0, 1.0]
blue =  [0.0, 0.0, 1.0, 1.0]

fp = open('EPOT.dat', 'w')

# Initialize alloc all variables

swobj = MDobj(strain, dirname, cohesv)
swobj.initialize()

nucleus = swobj.choose_elipse_state(np.array([0,0,0.38475]), 0.2, 0.2)

swobj.make_frk_dislocation(nucleus)
energy = swobj.eval(nucleus)
swobj.sw.writeatomeyecfg(dirname + '/img_B.cfg')
swobj.saveNucleus(nucleus, dirname + '/nucleus_B.cfg')
fp.write('Energy of state B = {0}\n'.format(energy))

print("Program exits normally")


