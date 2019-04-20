#from OpenGL.GLUT import *
#from OpenGL.GLU import *
#from OpenGL.GL import *
import sys
import numpy as np

sys.path.append('../../lib/') #axiom/lib contains mdsw
import lorenz
lrz = lorenz.Sim()

# Test run for continuing simulations
if 0:
    cp = [0,1,0]
    for i in range(100):
        print('epoch {0}'.format(i))
        lrz.run_fine_ODE(10, 10, 8./3, 25, 0, cp[0], cp[1], cp[2])
        traj = lrz.trajectory()
        cp = traj[-1,:]
        print('traj = ')
        print(traj)

cp = [0,1,0]

lrz.run_fine_ODE(10, 10, 8./3, 25, 0, cp[0], cp[1], cp[2])

traj = lrz.trajectory()




