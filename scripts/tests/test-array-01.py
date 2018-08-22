from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys

sys.path.append('/Users/x/Downloads/GraphicsScratch/axiom/lib/')
import mdsw
from View import Viewer

""" Main Program Starts
"""
sw = mdsw.SWFrame()
sw.initvars()
sw.dirname = "runs/si-example"
sw.crystalstructure = "diamond-cubic"
sw.latticeconst = mdsw.VectorDouble([5.4309529817532409, 5.4309529817532409, 5.4309529817532409])
sw.latticesize = mdsw.VectorDouble([ 1 ,0, 0, 4,  0, 1, 0, 4,  0, 0, 1, 4 ])
sw.makecrystal()
sw.finalcnfile = "relaxed.cn"
sw.writecn(0,False)
sw.eval()

y = sw.SR()
print(y)


print("Test fixed()")
f = sw.fixed()
print(f[0])

f[0] = 1

print(f[0])




