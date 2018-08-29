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

swobj = MDobj(strain, dirname)
swobj.reset()

# A set of test functions for MDobj class


#alg = GreedySaerch(swobj, cohesv, strain, dirname)

# A set of test functions for GreedySearch algorithm


#alg = RlSaerch(swobj, cohesv, strain, dirname)

# A set of test functions for Reinforcement Learning algorithm


# A search begins. 

#alg.search()

#alg.visualize_path()
  

