from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys

sys.path.append('/Users/x/Downloads/GraphicsScratch/axiom/lib/')
sys.path.append('/Users/x/Downloads/GraphicsScratch/axiom/scripts/tests/')
import mdsw
from View import Viewer

""" Conjugate-Gradient relaxation """
def relax_fixbox(obj,relaxsteps):
    obj.conj_ftol = 1e-7
    obj.conj_itmax = 1000
    obj.conj_fevalmax = relaxsteps
    obj.conj_fixbox = 1
    obj.relax()

""" Main Program Starts
"""
sw = mdsw.SWFrame()
sw.initvars()
sw.NNM = 200

epsilon = 0.052

nx = 12
ny = 12
nz = 16
sw.dirname = "runs/si-example"
sw.crystalstructure = "diamond-cubic"
sw.latticeconst = mdsw.VectorDouble([5.4309529817532409, 5.4309529817532409, 5.4309529817532409])
sw.latticesize = mdsw.VectorDouble([ 1 ,1, 0, nx,  -1, 1, 0, ny,  0, 0, 1, nz ])
sw.makecrystal()
sw.finalcnfile = "perf.cn"
sw.writecn(0,False)
sw.eval()


sw.saveH()
relax_fixbox(sw,sw.conj_fevalmax)

sw.vacuumratio = 1. - 1. /(1.+0.2)
sw.input = mdsw.VectorDouble([3,3,0.2])
sw.changeH_keepR()
sw.relax()
sw.finalcnfile = "0K_surf001.cn"
sw.writecn(0,False)

ny = sw.latticesize[7]
for i in range(int(ny)):
   ymin = -0.5006+1.0/ny*i
   ymax = -0.5006+1.0/ny*(i+0.5)
   sw.input = mdsw.VectorDouble([ 1, -10 , 10, ymin, ymax, 0.403,   10 ])
   sw.fixatoms_by_position()
   sw.input = mdsw.VectorDouble([ 1, -10 , 10, ymin, ymax  , -10 , -0.416 ])
   sw.fixatoms_by_position()

sw.input = mdsw.VectorDouble([1])
sw.setfixedatomsgroup()
sw.freeallatoms()

for i in range(int(ny)):
   ymin = -0.5006+1.0/ny*(i+0.5)
   ymax = -0.5006+1.0/ny*(i+1)
   sw.input = mdsw.VectorDouble([ 1, -10, 10 ,ymin, ymax, 0.403  , 10 ])
   sw.fixatoms_by_position()
   sw.input = mdsw.VectorDouble([ 1, -10, 10, ymin, ymax, -10 ,-0.416 ])
   sw.fixatoms_by_position()

sw.input = mdsw.VectorDouble([2])
sw.setfixedatomsgroup()
sw.freeallatoms()

sw.input = mdsw.VectorDouble([ 1,  0,  0.8, 0,  1 ])
sw.movegroup()
sw.input = mdsw.VectorDouble([ 1,  0, -0.8, 0,  2 ])
sw.movegroup()

relax_fixbox(sw, sw.conj_fevalmax)
sw.finalcnfile = "0K_0.0_relaxed_surf001.cn"
sw.writecn(0,False)
sw.finalcnfile = "0K_0.0_relaxed_surf001.cfg"
sw.writeatomeyecfg(sw.finalcnfile)

H = sw.H()
H11_0 = H[0][0]; H22_0 = H[1][1]; H33_0 = H[2][2]
sw.saveH()
H11_fix = H11_0*(1.0-epsilon)
H[0][0] = H11_fix #by reference, therefore sw.H has the same memory address as H
sw.conj_ftol = 1e-4
sw.conj_itmax = 3800
sw.conj_fevalmax = 6000
sw.conj_fixbox = 1
sw.relax()
sw.SHtoR()
sw.eval()
#set fp [ open "EPOT_1.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp
sw.finalcnfile = "0K_"+str(epsilon)+"_strained.cn"
sw.writecn(0,False)
sw.finalcnfile = "0K_"+str(epsilon)+"_strained.cfg"
sw.writeatomeyecfg(sw.finalcnfile)

NP = sw.NP
sw.refreshnnlist()
sw.input = mdsw.VectorDouble([NP])
sw.fprintnnlist()

view = Viewer(sw, 600, 600)
view.rendering()
#sw.sleep()



