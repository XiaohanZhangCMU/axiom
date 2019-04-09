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


#include "vector.hpp"

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
    Sim() { }

    int run();
};

void fillVBO();
void drawDirect();
void drawFigure();
void drawComponents();
void draw();
void timer(int valor);
void keyboardDown(unsigned char key, int x, int y);
void mouseClick(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void reshape(int w, int h);
void initGL(int width, int height);

} /* namespace axiom */
#endif /* __FEM1D_HPP__ */
