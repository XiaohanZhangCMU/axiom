from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

import numpy as np
import sys

class Viewer(object):
    def __init__(self, sim, width, height, display=None):
        #-----------
        # VARIABLES
        #-----------

        self.g_fViewDistance = 9.
        self.g_Width = 600   # Fix this to take args
        self.g_Height = 600  # Fix this to take args

        self.g_nearPlane = 1.
        self.g_farPlane = 1000.

        self.action = ""
        self.xStart = self.yStart = 0.
        self.zoom = 65.

        self.xRotate = 0.
        self.yRotate = 0.
        self.zRotate = 0.

        self.xTrans = 0.
        self.yTrans = 0.

        self.SR = sim.SR()
        self.fixed = sim.fixed()

        # For test only: plot only the first 100 fixed atoms
        self.fixed[:100] = 1
        self.SR = self.SR[self.fixed==1,:]
        self.atomradius = 0.05


    #-------------------
    # SCENE CONSTRUCTOR
    #-------------------

    def scenemodel(self):
        glRotate(90,0.,0.,1.)
        glutSolidTeapot(1.)

    def sphere(self):
        glPushMatrix(self)
        color = [1.0,0.,0.,1.]
        glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
        glutSolidSphere(2,20,20)
        glPopMatrix()

    def spheres_bruteforce(self):
        glPushMatrix()
        color = [1.0,0.,0.,1.]
        glTranslated(-1.2,-1.2,-1.2);
        glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
        glutSolidSphere(0.2,20,20)
        glPopMatrix()

        glPushMatrix()
        color = [1.0,0.,0.,1.]
        glTranslated(1.2,1.2,1.2);
        glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
        glutSolidSphere(0.2,20,20)
        glPopMatrix()


    def spheres(self):
        color = [1.0,0.,0.,1.]
        glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
        for ind in range(self.SR.shape[0]):
            glPushMatrix()
            glTranslated(self.SR[ind,0], self.SR[ind,1], self.SR[ind,2]);
            glutSolidSphere(self.atomradius,20,20)
            glPopMatrix()

    #--------
    # VIEWER
    #--------

    def printHelp(self):
        print("""\n\n
             -------------------------------------------------------------------\n
             Left Mousebutton       - move eye position (+ Shift for third axis)\n
             Middle Mousebutton     - translate the scene\n
             Right Mousebutton      - move up / down to zoom in / out\n
             key (r)                - reset viewpoint\n
             key (z)                - zoom in (not working yet)\n
             Key (q)                - exit the viewer\n
             -------------------------------------------------------------------\n
             \n""")


    def init(self):
        glEnable(GL_NORMALIZE)
        glLightfv(GL_LIGHT0,GL_POSITION,[ .0, 10.0, 10., 0. ] )
        glLightfv(GL_LIGHT0,GL_AMBIENT,[ .0, .0, .0, 1.0 ]);
        glLightfv(GL_LIGHT0,GL_DIFFUSE,[ 1.0, 1.0, 1.0, 1.0 ]);
        glLightfv(GL_LIGHT0,GL_SPECULAR,[ 1.0, 1.0, 1.0, 1.0 ]);
        glEnable(GL_LIGHT0)
        glEnable(GL_LIGHTING)
        glEnable(GL_DEPTH_TEST)
        glDepthFunc(GL_LESS)
        glShadeModel(GL_SMOOTH)
        self.resetView()


    def resetView(self):
        #global zoom, xRotate, yRotate, zRotate, xTrans, yTrans
        self.zoom = 65.
        self.xRotate = 0.
        self.yRotate = 0.
        self.zRotate = 0.
        self.xTrans = 0.
        self.yTrans = 0.
        glutPostRedisplay()


    def display(self):
        # Clear frame buffer and depth buffer
        glClearColor(1.0, 1.0, 1.0, 0.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        # Set up viewing transformation, looking down -Z axis
        glLoadIdentity()
        gluLookAt(0, 0, -self.g_fViewDistance, 0, 0, 0, -.1, 0, 0)   #-.1,0,0
        # Set perspective (also zoom)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(self.zoom, float(self.g_Width)/float(self.g_Height), self.g_nearPlane, self.g_farPlane)
        glMatrixMode(GL_MODELVIEW)
        # Render the scene
        self.polarView()

        #scenemodel()
        #sphere()
        self.spheres()

        # Make sure changes appear onscreen
        glutSwapBuffers()


    def reshape(self, width, height):
        #global g_Width, g_Height
        self.g_Width = width
        self.g_Height = height
        glViewport(0, 0, self.g_Width, self.g_Height)


    def polarView(self):
        glTranslatef( self.yTrans/100., 0.0, 0.0 )
        glTranslatef(  0.0, -self.xTrans/100., 0.0)
        glRotatef( -self.zRotate, 0.0, 0.0, 1.0)
        glRotatef( -self.xRotate, 1.0, 0.0, 0.0)
        glRotatef( -self.yRotate, .0, 1.0, 0.0)


    def keyboard(self, key, x, y):
        global zTr, yTr, xTr
        if(key==b'r'): self.resetView()
        if(key==b'z'): self.zoomin()
        if(key==b'q'): sys.exit(0)
        glutPostRedisplay()


    def mouse(self, button, state, x, y):
        #global action, xStart, yStart
        if (button==GLUT_LEFT_BUTTON):
            if (glutGetModifiers() == GLUT_ACTIVE_SHIFT):
                self.action = "MOVE_EYE_2"
            else:
                self.action = "MOVE_EYE"
        elif (button==GLUT_MIDDLE_BUTTON):
            self.action = "TRANS"
        elif (button==GLUT_RIGHT_BUTTON):
            self.action = "ZOOM"
        self.xStart = x
        self.yStart = y


    def motion(self, x, y):
        #global zoom, xStart, yStart, xRotate, yRotate, zRotate, xTrans, yTrans
        if (self.action=="MOVE_EYE"):
            self.xRotate += x - self.xStart
            self.yRotate -= y - self.yStart
        elif (self.action=="MOVE_EYE_2"):
            self.zRotate += y - self.yStart
        elif (self.action=="TRANS"):
            self.xTrans += x - self.xStart
            self.yTrans += y - self.yStart
        elif (self.action=="ZOOM"):
            self.zoom -= y - self.yStart
            if self.zoom > 150.:
                self.zoom = 150.
            elif self.zoom < 1.1:
                self.zoom = 1.1
        else:
            print("unknown action\n", action)
        self.xStart = x
        self.yStart = y
        glutPostRedisplay()

    def zoomin(self):
        scale_fac = 1.1
        glPushMatrix()
        glScalef(scale_fac,scale_fac,scale_fac)
        glPopMatrix()
        #self.zoom -= 0 - self.yStart
        #if self.zoom > 150.:
        #    self.zoom = 150.
        #elif self.zoom < 1.1:
        #    self.zoom = 1.1
        glutPostRedisplay()


    #------
    # MAIN
    #------
    def rendering(self):
        # GLUT Window Initialization
        glutInit()
        glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB| GLUT_DEPTH)      # zBuffer
        glutInitWindowSize (self.g_Width,self.g_Height)
        glutInitWindowPosition (0 + 4, int(self.g_Height / 4) )
        glutCreateWindow ("Visualizzatore_2.0")
        # Initialize OpenGL graphics state
        self.init ()
        # Register callbacks
        glutReshapeFunc(self.reshape)
        glutDisplayFunc(self.display)
        glutMouseFunc(self.mouse)
        glutMotionFunc(self.motion)
        glutKeyboardFunc(self.keyboard)
        self.printHelp()
        # Turn the flow of control over to GLUT
        glutMainLoop()
