/*
 *  solver.h
 */
#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <vector>
#include <string>
#include "linalg3.h"

using namespace std;

class Lorenz {
public:
	Lorenz(Vector3 txyz0) { this->txyz = this->txyz0=txyz0; }
	virtual Vector3 f(double t, Vector3 txyz) { return txyz0 ; } 
    virtual ~Lorenz() {} 

	Vector3 txyz0; //initial phase point 
	Vector3 txyz;  //current phase point
	double x, y, z;
};

class fLorenz : public Lorenz {
public:
	fLorenz(double sigma, double r, double b, Vector3 xyz0) : Lorenz(xyz0) 
	{ this->sigma = sigma; this->b = b; this->r = r; }
	Vector3 f(double t, Vector3 txyz);
	double sigma, b, r;
};

class fLorenz84 : public Lorenz {
public:
	fLorenz84(double sigma, double r, double b, Vector3 xyz0) : Lorenz(xyz0) 
	{ this->sigma = sigma; this->b = b; this->r = r; }
	Vector3 f(double t, Vector3 txyz);
	double sigma, b, r;
};

class fPickover : public Lorenz {
public:
	fPickover(double sigma, double r, double b, Vector3 xyz0) : Lorenz(xyz0) 
	{ this->sigma = sigma; this->b = b; this->r = r; }
	Vector3 f(double t, Vector3 txyz);
	double sigma, b, r;
};

class fRabinovich : public Lorenz {
public:
	fRabinovich(double gam, double alphavalue, Vector3 xyz0) : Lorenz(xyz0) 
	{ this->gam = gam; this->alphavalue = alphavalue; }
	Vector3 f(double t, Vector3 txyz);
	double gam, alphavalue;
};

class fRossler : public Lorenz {
public:
	fRossler(double a, double b, double c, Vector3 xyz0) : Lorenz(xyz0) 
	{ this->a = a; this->b = b; this->c = c; }
	Vector3 f(double t, Vector3 txyz);
	double a, b, c;
};

class Integrator 
{
public:
	Integrator(Lorenz* lrz) { this->lrz = lrz; }
	virtual void step(double t, Vector3 &y, double h) { }
    virtual ~Integrator() { if (!lrz) delete lrz; } 
	inline Vector3 f(double t, Vector3 xyz) { return  this->lrz->f(t,xyz); }  
	Lorenz *lrz;
};
class IntegratorEuler : public Integrator {
public:
	IntegratorEuler(Lorenz* lrz) : Integrator(lrz) { } ;
	void step(double t, Vector3 &y, double h);	
};
class IntegratorEulerMod : public Integrator {
public:
	IntegratorEulerMod(Lorenz* lrz) : Integrator(lrz) { } ;
	void step(double t, Vector3 &y, double h);	
};
class IntegratorRungeKutta2 : public Integrator {
public:
	IntegratorRungeKutta2(Lorenz* lrz) : Integrator(lrz) { } ;
	void step(double t, Vector3 &y, double h);	
};
class IntegratorRungeKutta3 : public Integrator {
public:
	IntegratorRungeKutta3(Lorenz* lrz) : Integrator(lrz) { } ;
	void step(double t, Vector3 &y, double h);	
};
class IntegratorRungeKutta4 : public Integrator {
public:
	IntegratorRungeKutta4(Lorenz* lrz) : Integrator(lrz) { } ;
	void step(double t, Vector3 &y, double h);	
};
class IntegratorRungeKutta45 : public Integrator {
public:
	IntegratorRungeKutta45(Lorenz* lrz) : Integrator(lrz) { } ;
	void step(double t, Vector3 &y, double h);	
};

#define min(a,b) ((a<b)?(a):(b))
#define max(a,b) ((a>b)?(a):(b))
#define clamp(mn,v,mx) max(min(v,mx),mn)

#endif /* __SOLVER_H__ */
//void step(double &t, Vector3 &y, double &h);
