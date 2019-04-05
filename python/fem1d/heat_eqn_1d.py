#from OpenGL.GLUT import *
#from OpenGL.GLU import *
#from OpenGL.GL import *
import sys
import numpy as np

sys.path.append('../../lib/') #axiom/lib contains mdsw
import fem

def exact (  x ):
  return x * ( 1.0 - x ) * np.exp ( x );

n = 11
a = 0.0
b = 1.0
ua = 0.0
ub = 0.0
ob = fem.Fem()
ob.set_mesh( n, a, b )
ob.fem1d_heat_steady ( n, a, b, ua, ub)
u = ob.getU()
X = ob.getX()

print("FEM solution = ")
print(u)

exact_solution = [ ]
for xx in X:
    exact_solution.append(exact(xx))

print("Exact solution = ")
print(np.array(exact_solution))
