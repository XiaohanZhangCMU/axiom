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
#from SearchAlgorithms import GreedySearch
from SearchAlgorithms import DQNSearch

sys.path.append('../../lib/')
sys.path.append('../../tests/')
import mdsw
#from View import Viewer

""" Specify a small B state, check if DQN works as expected.
    Run this script. Check EPOT.dat and path images.
"""

# Control parameters

strain = "0.05"
dirname = 'dqn_neighbor_search_'+ strain +'/'

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

# Make sure loop through the plane having more atoms around stateB
# Or loop through the bottom plane. Loop through up plane has a problem:
# the last two atoms cannot be found at the end of an episode
# This needs to be fixed !!!!!!
stateA_u = np.intersect1d(swobj.pairs[:,0], stateA)
stateB_u = np.intersect1d(swobj.pairs[:,1], stateB)
stateA_d = np.intersect1d(swobj.pairs[:,0], stateA)
stateB_d = np.intersect1d(swobj.pairs[:,1], stateB)
print(len(stateA_u))
print(len(stateA_d))
print(len(stateB_u))
print(len(stateB_d))
stateB = stateB_u if len(stateB_u)>len(stateB_d) else stateB_d
stateA = stateA_u if len(stateA_u)>len(stateA_d) else stateA_d

alg = DQNSearch(stateA, stateB,
                learning_rate=0.1,
                reward_decay=0.9,
                e_greedy=0.9,
                replace_target_iter=20,
                memory_size=2000,
                output_graph=True
               )

# Begin search
alg.train(swobj, n_episodes = 300)
alg.write_path(swobj)


'''
env = Maze()
alg = DQNSearch(cohesv, strain, dirname,
                env.n_actions, env.n_features,
                learning_rate=0.01,
                reward_decay=0.9,
                e_greedy=0.9,
                replace_target_iter=200,
                memory_size=2000,
                output_graph=True
               )
alg.search(300, env)
'''


print("Program exits normally")


