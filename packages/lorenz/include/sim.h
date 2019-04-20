#ifndef __SIM_HPP__
#define __SIM_HPP__

#include <fcntl.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <OpenGL/OpenGL.h>
#include <OpenGL/gl.h>
#include <OpenGL/CGLDevice.h>
#include <GLUT/glut.h>


#include "linalg3.h"

#include "drawing.h"
#include "solver.h"
#include "vbo.h"

#define TIMER_DELAY 25
#define WIDTH 1280
#define HEIGHT 800
#define POSX 10
#define POSY 10

#define MAX_POINTS 1000000
#define MAX_SKIP 10

#define NUM_SHOW 1000


namespace axiom {

static int lastx=0;
static int lasty=0;

class Sim {
    public:
// Camera
    Sim() : trajectory(0) {
      //  this->max_points = MAX_POINTS;
        this->dt = 0.001;
        this->t = 0;
        this->points = MAX_POINTS;
    }

    ~Sim() { 
//        if (trajectory) free(this->trajectory);
    }

    void run_fine_ODE(int points, double sigma, double beta,
            double gamma, double t0, double x0, double y0, double z0);

    class Vector3 *trajectory;
    int points;
    double dt;
    double t;
};

//void drawDirect();
//void drawFigure();
//void drawComponents();
//void draw();
//void timer(int valor);
//void keyboardDown(unsigned char key, int x, int y);
//void mouseClick(int button, int state, int x, int y);
//void mouseMotion(int x, int y);
//void reshape(int w, int h);
//void initGL(int width, int height);

} /* namespace axiom */
#endif /* __FEM1D_HPP__ */
