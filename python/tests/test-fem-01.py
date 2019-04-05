import sys
sys.path.append('../..//lib/')
import mdfem


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


fem = mdfem.FEMFrame()
fem.contact_file = "new"
v = mdfem.VectorInt([1,1])
fem.bNds = v
fem.bTags = v
fem.bNds[0] = 0
fem.incnfile = "newfile"
fem.finalcnfile = "oldfile"
n = 10
fem.setnolog(n)
md = mdfem.MDFrame()
md.eval()


