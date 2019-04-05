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


""" Create state A, A0 and a series of B states and evaluate.
    Run this script. Check EPOT.dat and img_Bx.cfg files.
"""

# Control parameters

strain = "0.05"
dirname = '/home/x/runs/test_state_AnB_energy_'+ strain +'/'

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


# Initialize alloc all variables and create State A0 under strain:

swobj = MDobj(strain, dirname, cohesv)
swobj.initialize()

fp = open(swobj.dirname+'/EPOT.dat', 'w')
fp.write('Energy of state A0 = {0}\n'.format(swobj.E0))

#H = swobj.sw.H
#H[0] = H[0]*(1.0-np.double(strain))
#swobj.sw.H = H
#swobj.relax_fixbox()
#swobj.sw.SHtoR()
#energy = swobj.eval( np.array([]))
#swobj.sw.finalcnfile = dirname+"0K_0.0_relaxed_surf001.cn"
#swobj.sw.writecn(0,False)
#swobj.sw.finalcnfile = dirname+"0K_0.0_relaxed_surf001.cfg"
#swobj.sw.writeatomeyecfg(swobj.sw.finalcnfile)
#swobj.restoreConfig()

# Create State A:

nucleus = swobj.choose_elipse_state(np.array([0,0,0.38475]), 0.08, 0.08)

swobj.make_frk_dislocation(nucleus)
energy = swobj.eval(nucleus)
swobj.sw.writeatomeyecfg(dirname + '/img_A.cfg')
swobj.saveNucleus(nucleus, dirname + '/nucleus_A.cfg')
fp.write('Energy of state A = {0}\n'.format(energy))

swobj.restoreConfig()

# Create a series of State B1:

Bcount = 0
for radius in np.linspace(0.1, 0.4, 7):
    nucleus = swobj.choose_elipse_state(np.array([0,0,0.38475]), radius, radius)
    swobj.make_frk_dislocation(nucleus)
    energy = swobj.eval(nucleus)
    swobj.sw.writeatomeyecfg(dirname + '/img_B'+str(Bcount)+'.cfg')
    swobj.saveNucleus(nucleus, dirname + '/nucleus_B'+str(Bcount)+'.cfg')
    fp.write('Energy of state B'+str(Bcount)+' = {0}\n'.format(energy))
    swobj.restoreConfig()
    Bcount += 1

swobj.restoreConfig()
fp.close()

print("Program exits normally")


