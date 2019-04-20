# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

#include "sim.h"

namespace axiom {

//double rotx = 0;
//double roty = 0;
//double zoom = 0.8;
//double autorotate = 0.2;
//int screenWidth = WIDTH;
//int screenHeight = HEIGHT;

//int fps = 0;
//char fps_str[10]="";
//char param_str[500]="";

// Status
//bool running = 1;
//int showAll = 1;
//int showComponents = 1;
//int lineMode = 1;
//int moveIncrement = 100;

//int counter = 0;

//int usingVBO=1;

//double minx, maxx, miny, maxy, minz, maxz;
//double alpha=0.8;
//int skip=0;
//int points = 1000;
//int solver;
//Vector currentPoint(0,1,-2,1);
//Vector initial(0,1,-2,1);
//double SIGMA = 10;
//double B = 8.0/3.0;
//double R = 28;

    #define Realloc(p,t,s) {p=(t *)realloc(p,sizeof(t)*s);}

void Sim::run_fine_ODE(int points, double sigma, double beta, double gamma, double t0, double x0, double y0, double z0) {
	
	printf("Filling VBO\n");

    this->points = points;
    //if (trajectory != NULL)
    //    free(trajectory);
    trajectory = (Vector3 *) realloc(trajectory, points*sizeof(Vector3));
	//printf("I am here 1\n");

//	COLOUR c;
//
    Vector3 currentPoint(x0, y0, z0);
//
	Lorenz * lrz = new fLorenz(sigma, gamma, beta, currentPoint);
	Integrator *integrator = new IntegratorEuler(lrz);

	//printf("I am here 2\n");

//	minx = miny = minz = 10000000000;
//	maxx = maxy = maxz = -10000000000;
//	
	for (int i = 0; i < points; i++) {
		integrator->step(t, currentPoint, dt);
	    //printf("I am here 2.1\n");

        trajectory[i] = currentPoint;
	    //printf("I am here 2.2\n");
        //currentPoint[0] += dt;
	    //printf("I am here 2.3\n");
		
		//vertices[i].position[0] = currentPoint[1];
		//vertices[i].position[1] = currentPoint[2];
		//vertices[i].position[2] = currentPoint[3];
		//
		//minx = min(vertices[i].position[0], minx);
		//maxx = max(vertices[i].position[0], maxx);

		//miny = min(vertices[i].position[1], miny);
		//maxy = max(vertices[i].position[1], maxy);
		//
		//minz = min(vertices[i].position[2], minz);
		//maxz = max(vertices[i].position[2], maxz);

		std::cout<<"cp = ("<< currentPoint.x<<", "<<currentPoint.y<<", "<<currentPoint.z<<")"<<std::endl;
	}
//	
//	updateVBO(vertices, points);
//	
    if (!lrz) delete(lrz);
	if (!integrator) delete(integrator);
	printf("VBO Filled OK\n");
}
}
