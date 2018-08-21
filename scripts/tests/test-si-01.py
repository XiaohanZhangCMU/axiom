import sys
sys.path.append('/Users/x/Downloads/GraphicsScratch/axiom/lib/')
import mdsw

def relax_fixbox(obj,relaxsteps):
    #Conjugate-Gradient relaxation
    obj.conj_ftol = 1e-7
    obj.conj_itmax = 1000
    obj.conj_fevalmax = relaxsteps
    obj.conj_fixbox = 1
    obj.relax()

def open_window(obj):
    obj.atomradius = mdsw.VectorDouble([0.67])
    obj.bondradius = 0.3
    obj.bondlength = 2.8285 #for Si
    obj.atomcolor = [ "orange"]
    obj.highlightcolor = "purple"
    obj.bondcolor = "red"
    obj.backgroundcolor = "gray70"
    obj.plotfreq = 1
    obj.rotateangles = mdsw.VectorDouble([0,0,0,1.25])
    obj.plot_limits = mdsw.VectorDouble([  1, -10 , 10 , -10 , 10 , 0.15 , 0.5  ])
    obj.refreshnnlist
    obj.openwin()
    obj.alloccolors
    obj.rotate
    obj.saverot
    obj.plot_atom_info = 0
    obj.plot_color_windows=mdsw.VectorDouble([1, -10, 10])
    obj.plot()

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

open_window(sw)
sw.sleep()
relax_fixbox(sw,1000)



