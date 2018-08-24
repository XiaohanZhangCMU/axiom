from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys

sys.path.append('../../lib/')
import mdsw
from View import Viewer

""" Conjugate-Gradient relaxation """
def relax_fixbox(obj,relaxsteps):
    obj.conj_ftol = 1e-7
    obj.conj_itmax = 1000
    obj.conj_fevalmax = relaxsteps
    obj.conj_fixbox = 1
    obj.relax()

#def open_window(obj):
#    obj.atomradius = mdsw.VectorDouble([0.67])
#    obj.bondradius = 0.3
#    obj.bondlength = 2.8285 #for Si
#    obj.atomcolor = [ "orange"]
#    obj.highlightcolor = "purple"
#    obj.bondcolor = "red"
#    obj.backgroundcolor = "gray70"
#    obj.plotfreq = 1
#    obj.rotateangles = mdsw.VectorDouble([0,0,0,1.25])
#    obj.plot_limits = mdsw.VectorDouble([  1, -10 , 10 , -10 , 10 , 0.15 , 0.5  ])
#    obj.refreshnnlist
#    obj.openwin()
#    obj.alloccolors
#    obj.rotate
#    obj.saverot
#    obj.plot_atom_info = 0
#    obj.plot_color_windows=mdsw.VectorDouble([1, -10, 10])
#    obj.plot()

#def render(name):
#    glutInit(sys.argv)
#    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
#    glutInitWindowSize(400,400)
#    glutCreateWindow(name)
#
#    glClearColor(0.,0.,0.,1.)
#    glShadeModel(GL_SMOOTH)
#    glEnable(GL_CULL_FACE)
#    glEnable(GL_DEPTH_TEST)
#    glEnable(GL_LIGHTING)
#    lightZeroPosition = [10.,4.,10.,1.]
#    lightZeroColor = [0.8,1.0,0.8,1.0] #green tinged
#    glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
#    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
#    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1)
#    glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05)
#    glEnable(GL_LIGHT0)
#    glutDisplayFunc(display)
#    glMatrixMode(GL_PROJECTION)
#    gluPerspective(40.,1.,1.,40.)
#    glMatrixMode(GL_MODELVIEW)
#    gluLookAt(0,0,10,
#              0,0,0,
#              0,1,0)
#    glPushMatrix()
#    glutMainLoop()
#    return
#
#def display():
#    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
#    glPushMatrix()
#    color = [1.0,0.,0.,1.]
#    glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
#    glutSolidSphere(2,20,20)
#    glPopMatrix()
#    glutSwapBuffers()
#    return


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

#render("mdsw")
#open_window(sw)
view = Viewer(sw, 600, 600)
view.rendering()
sw.sleep()
relax_fixbox(sw,1000)



