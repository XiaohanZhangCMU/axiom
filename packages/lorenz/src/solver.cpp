/*
 *  lorenz
 */

#include <math.h>
#include "linalg3.h"
#include "solver.h"

//char *solverNames[] = { "Euler", "EulerMod", "RungeKutta2", "RungeKutta3", "RungeKutta4", "RKF" };

//int skip = 2;
//int points = 100000;

//double sigma; //= 10;
//double b; //= 8.0 / 3.0;
//double r; //= 28;

Vector3 fLorenz::f(double t, Vector3 txyz) {
    	Vector3 f;
    	x = txyz[0]; y=txyz[1]; z=txyz[2];
    	f[0] = - sigma * x + sigma * y;
    	f[1] = - x * z + r * x - y;
    	f[2] = x * y - b * z;
//	std::cout<<"sigma ="<< sigma <<"; r = "<<r<<"; b = "<<b<<std::endl;
//	std::cout<<"f = ("<< f[0]<<", "<<f[1]<<", "<<f[2]<<", "<<f[3]<<")"<<std::endl;
//	std::cout<<"xyz = ("<< x<<", "<<y<<", "<<z<<")"<<std::endl;
    	return f;
}

Vector3 fLorenz84::f(double t, Vector3 txyz) {
	
	double A=sigma;
	double B=b;
	double F=r;
	double G=1;
	
    	Vector3 f;
    	x = txyz[0]; y=txyz[1]; z=txyz[2];
    	f[0] = - A * x - pow(y, 2) - pow(z, 2) + A*F;
    	f[1] = -y + x*y - B * x * z + G; 
    	f[2] = -z + B * x*y + x*z;
    	return f;
}

Vector3 fPickover::f(double t, Vector3 txyz) {
	
	sigma = 1;
	b = 1.8;
	r = 0.71;
	
	double A = sigma;
	double B = b;
	double C = r;
	double D = 1.51;
	
	Vector3 f;
    	x = txyz[0]; y=txyz[1]; z=txyz[2];
	f[0] = sin ( A*y) - z * cos(B * x);
	f[1] = z * sin (C*x) - cos(D * y);
	f[2] = sin(x);
	return f;
}

//double gam = 0.87;
//double alphavalue = 1.1;

// http://en.wikipedia.org/wiki/Rabinovich-Fabrikant_equations
Vector3 fRabinovich::f(double t, Vector3 txyz) {
        double r = sqrt(x*x + y*y);
    	Vector3 f;
    	x = txyz[0]; y=txyz[1]; z=txyz[2];
    	f[0] = y*(z-1+x*x)+gam*x;
    	f[1] = x*(3*z+1-x*x)+gam*y;
    	f[2] = -2*z*(alphavalue+x*y);
    	return f;
}


// http://en.wikipedia.org/wiki/R%C3%B6ssler_map
Vector3 fRossler::f(double t, Vector3 txyz) {
	    double a = 0.1;
	    double b = 0.1;
	    double c = 14;
	    
	    double r = sqrt(x*x + y*y);
	
    	Vector3 f;
    	x = txyz[0]; y=txyz[1]; z=txyz[2];
    	f[0] = -y-z;
    	f[1] = x+a*y;
    	f[2] = b+z*(x-c);
    	return f;
}

void IntegratorEuler::step(double t, Vector3 &y , double h) {
//	std::cout<<"P:y[0] = "<<y[0]<<std::endl;
	y+= h*f(t, y);
//	std::cout<<"A:y[0] = "<<y[0]<<std::endl;
}

void IntegratorEulerMod::step(double t, Vector3 &y , double h) {
	Vector3 k1 = h * f(t, y);
	Vector3 k2 = h * f(t + 0.5*h, y + 0.5 * k1);
	y+= k2;
}

void IntegratorRungeKutta2::step(double t, Vector3 & y, double h) {
	Vector3 k1 = h * f(t, y);
	Vector3 k2 = h * f(t + 0.5*h, y + 0.5 * k1);
	
	y+= k2;
}

void IntegratorRungeKutta3::step(double t, Vector3& y, double h) {
    	Vector3 k1 = h * f(t, y);
    	Vector3 k2 = h * f(t + 0.5*h, y + 0.5 * k1);
    	Vector3 k3 = h * f(t + h, y -k1 + 2* k2);
    	    
    	y += (k1 + 4*k2 + k3) / 6.0;
}


void IntegratorRungeKutta4::step(double t, Vector3& y, double h) {
    	Vector3 k1 = h * f(t, y);
    	Vector3 k2 = h * f(t + 0.5*h, y + 0.5 * k1);
    	Vector3 k3 = h * f(t + 0.5*h, y + 0.5 * k2);
    	Vector3 k4 = h * f(t+h, y + k3);
    	y += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

void IntegratorRungeKutta45::step(double t, Vector3& y, double h) {
// Runge Kutta
double error = 10000;
double minh = 0.0001;
double maxh = 0.01;
double dt = minh;

    	Vector3 k1 = h * f(t, y);
    	Vector3 k2 = h * f(t + 0.25*h, y + 0.25 * k1);
    	Vector3 k3 = h * f(t + (3.0/8.0)*h, y + (3.0/32.0)*k1 + (9.0/32.0) * k2);
    	Vector3 k4 = h * f(t + (12.0/13.0)*h, y + (1932.0/2197.0)*k1 - (7200.0/2197.0)*k2 + (7296.0/2197.0) * k3);
	Vector3 k5 = h * f(t + h, y + (439.0/216.0)*k1 - 8*k2 + (3680.0/513.0) * k3 - (845.0/4104.0) * k4);
	Vector3 k6 = h * f(t + 0.5*h, y - (8.0/27.0)*k1 + 2*k2 - (3544.0/2565.0)*k3 + (1859.0/4104.0)*k4 - (11.0/40.0)*k5);
	
	Vector3 ord4;
	Vector3 ord5;
	ord4 = y + (25.0/216.0)*k1 + (1408.0/2565.0)*k3 + (2197/4101)*k4 - (1.0/5.0)*k5;
	ord5 = y + (16.0/135.0)*k1 + (6656.0/12825.0)*k3 + (28561.0/56430.0)*k4 - (9.0/50.0)*k5 + (2.0/55.0) * k6;
	
	y = ord4;
	
	double s = pow( (error * h * 0.5) / (ord5-ord4).norm(), 0.25);
	
	h = clamp(minh, s * h, maxh);
	//std::cout << ord4 << " ¡ " << ord5 << " ¡ " << (ord5-ord4).abs() << " s="<< s << " h="<<h<<std::endl;
}
