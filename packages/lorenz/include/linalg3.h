/*
  linalg3.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sun Dec  3 16:06:53 2006

  FUNCTION  :  Linear Algebra package for 3 dimension array and matrix

  Featuring :  1. 3x3 Matrix
               2. 3 Vector
               3. UnitCell
*/

#ifndef _LINALG3_H
#define _LINALG3_H
#include <math.h>

class Vector3
{
public:
    double x, y, z;
    Vector3():x(0),y(0),z(0){};
    double operator [] (int i) const { return (&x)[i]; }
    double & operator [] (int i) { return (&x)[i]; }
    Vector3(double a,double b,double c):x(a),y(b),z(c){};
    Vector3(const double a[]){x=a[0];y=a[1];z=a[2];}
    
    //Vector3 operator + (const Vector3 &a) const
    //{ Vector3 c; c.x=x+a.x;c.y=y+a.y;c.z=z+a.z; return c; }
    //Vector3 operator - (Vector3 &a) const
    //{ Vector3 c; c.x=x-a.x;c.y=y-a.y;c.z=z-a.z; return c; }
    //Vector3 operator * (double b) const
    //{ Vector3 c; c.x=x*b;c.y=y*b;c.z=z*b; return c; }
    //Vector3 operator / (double b) const
    //{ Vector3 c; c.x=x/b;c.y=y/b;c.z=z/b; return c; }

    
    void addnv(double a,const Vector3& b){x+=a*b.x;y+=a*b.y;z+=a*b.z;}
    void invert(){x=1./x;y=1./y;z=1./z;}
    void subint() {x=x-rint(x);y=y-rint(y);z=z-rint(z);}
    double norm2() const {return x*x+y*y+z*z; }
    double norm() const  {return sqrt(x*x+y*y+z*z); }
    void clear() {x=y=z=0;}
    
    double operator * (const Vector3 &a) const
    { return (x*a.x+y*a.y+z*a.z); }
    
#if 0 /* return value */
    Vector3 & operator += (double b) { x+=b;y+=b;z+=b; return (*this); }
    Vector3 & operator -= (double b) { x-=b;y-=b;z-=b; return (*this); }
    Vector3 & operator *= (double b) { x*=b;y*=b;z*=b; return (*this); }
    Vector3 & operator /= (double b) { x/=b;y/=b;z/=b; return (*this); }
    Vector3 & operator += (const Vector3 &a) { x+=a.x;y+=a.y;z+=a.z; return (*this); }
    Vector3 & operator -= (const Vector3 &a) { x-=a.x;y-=a.y;z-=a.z; return (*this); }
    Vector3 & operator *= (const Vector3 &a) { x*=a.x;y*=a.y;z*=a.z; return (*this); }
    Vector3 & operator /= (const Vector3 &a) { x/=a.x;y/=a.y;z/=a.z; return (*this); }
#else /* do not return value */
    void operator += (double b) { x+=b;y+=b;z+=b; }
    void operator -= (double b) { x-=b;y-=b;z-=b; }
    void operator *= (double b) { x*=b;y*=b;z*=b; }
    void operator /= (double b) { x/=b;y/=b;z/=b; }
    void operator += (const Vector3 &a) { x+=a.x;y+=a.y;z+=a.z; }
    void operator -= (const Vector3 &a) { x-=a.x;y-=a.y;z-=a.z; }
    void operator *= (const Vector3 &a) { x*=a.x;y*=a.y;z*=a.z; }
    void operator /= (const Vector3 &a) { x/=a.x;y/=a.y;z/=a.z; }
#endif
    
    void set(const double a,const double b,const double c){x=a;y=b;z=c;}
    void set(const double a[]){x=a[0];y=a[1];z=a[2];}
    void add(const double a,const double b,const double c){x+=a;y+=b;z+=c;}
    void copytoarray(double a[]) const {a[0]=x;a[1]=y;a[2]=z;}
    Vector3 sq() const { Vector3 c; c.x=x*x;c.y=y*y;c.z=z*z; return c;}
    Vector3 sqroot() const { Vector3 c; c.x=(double)sqrt(x);c.y=sqrt(y);c.z=sqrt(z); return c;}
    
    Vector3 & orth(Vector3 &a)
    {
        double sub,norm2;
        norm2=a.norm2();
        if(norm2>1e-20)
        {
            sub=((*this)*a)/norm2;
            x-=a.x*sub;y-=a.y*sub;z-=a.z*sub;
        }
        return (*this);
    }
    Vector3 & proj(Vector3 &a)
    {
        double sub,norm2;
        norm2=a.norm2();
        if(norm2>1e-20)
        {
            sub=((*this)*a)/norm2;
            x=a.x*sub;y=a.y*sub;z=a.z*sub;
        }
        else clear();
        return (*this);
    }
    friend double dot(Vector3 &a, Vector3 &b)
    {
        double c;
        c=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        return c;
    }
    friend Vector3 cross(Vector3 &a, Vector3 &b)
    {
        Vector3 c;
        c[0]=a[1]*b[2]-a[2]*b[1];
        c[1]=a[2]*b[0]-a[0]*b[2];
        c[2]=a[0]*b[1]-a[1]*b[0];
        return c;
    }
};

extern Vector3 operator + (const Vector3 &a, const Vector3 &b);  
extern Vector3 operator - (const Vector3 &a, const Vector3 &b);
extern Vector3 operator * (const Vector3 &dv, const double d);
extern Vector3 operator * (const double d, const Vector3 &dv); 
extern Vector3 operator / (const Vector3 &dv, const double d); 

#endif // _LINALG3_H

