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
    
    Vector3 operator + (const Vector3 &a) const
    { Vector3 c; c.x=x+a.x;c.y=y+a.y;c.z=z+a.z; return c; }
    Vector3 operator - (Vector3 &a) const
    { Vector3 c; c.x=x-a.x;c.y=y-a.y;c.z=z-a.z; return c; }
    Vector3 operator * (double b) const
    { Vector3 c; c.x=x*b;c.y=y*b;c.z=z*b; return c; }
    Vector3 operator / (double b) const
    { Vector3 c; c.x=x/b;c.y=y/b;c.z=z/b; return c; }

    void add     (Vector3 &a, Vector3& b){x=a.x+b.x; y=a.y+b.y; z=a.z+b.z;}
    void subtract(Vector3 &a, Vector3& b){x=a.x-b.x; y=a.y-b.y; z=a.z-b.z;}
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
    friend Vector3 schmidt(Vector3 &a, Vector3 &b)
    {
#ifndef TINYNORM2
#define TINYNORM2 1e-10
#endif
        Vector3 c,a1;
        a1=a;
        c=a1;c.orth(b);
        if(c.norm2()<TINYNORM2) a1[0]+=1;
        c=a1;c.orth(b);
        if(c.norm2()<TINYNORM2) a1[1]+=1;
        c=a1;c.orth(b);
        if(c.norm2()<TINYNORM2)
        {
            FATAL("Schmidt("<<a<<", "<<b<<"): "
                  <<"too small! c="<<c);
        }
#undef TINYNORM
        return c/c.norm();
    }
    friend LOStream & operator <<(LOStream &os, const Vector3 &a)
    {//formatted output
        os.Printf("%18.10e %18.10e %18.10e",a.x,a.y,a.z);
        return os;
    }
    friend LIStream & operator >>(LIStream &is, Vector3 &a)
    {
        is >> a.x >> a.y >> a.z;
        return is;
    }

};

class Matrix33
{
public:
    double element[3][3];
public:
    Matrix33(){clear();}
    Matrix33(double e1,double e2,double e3,
             double e4,double e5,double e6,
             double e7,double e8,double e9)
        {element[0][0]=e1;element[0][1]=e2;element[0][2]=e3;
        element[1][0]=e4;element[1][1]=e5;element[1][2]=e6;
        element[2][0]=e7;element[2][1]=e8;element[2][2]=e9;}
    Matrix33(double e[9])
        {element[0][0]=e[0];element[0][1]=e[1];element[0][2]=e[2];
        element[1][0]=e[3];element[1][1]=e[4];element[1][2]=e[5];
        element[2][0]=e[6];element[2][1]=e[7];element[2][2]=e[8];}
    void set(const double e1,const double e2,const double e3,
             const double e4,const double e5,const double e6,
             const double e7,const double e8,const double e9)
        {element[0][0]=e1;element[0][1]=e2;element[0][2]=e3;
        element[1][0]=e4;element[1][1]=e5;element[1][2]=e6;
        element[2][0]=e7;element[2][1]=e8;element[2][2]=e9;}
    void set(const double e[])
        {element[0][0]=e[0];element[0][1]=e[1];element[0][2]=e[2];
        element[1][0]=e[3];element[1][1]=e[4];element[1][2]=e[5];
        element[2][0]=e[6];element[2][1]=e[7];element[2][2]=e[8];}
    void setcol(Vector3 &a,Vector3 &b,Vector3 &c)
        {element[0][0]=a[0];element[0][1]=b[0];element[0][2]=c[0];
        element[1][0]=a[1];element[1][1]=b[1];element[1][2]=c[1];
        element[2][0]=a[2];element[2][1]=b[2];element[2][2]=c[2];}
    void setcol(double a[],double b[],double c[])
        {element[0][0]=a[0];element[0][1]=b[0];element[0][2]=c[0];
        element[1][0]=a[1];element[1][1]=b[1];element[1][2]=c[1];
        element[2][0]=a[2];element[2][1]=b[2];element[2][2]=c[2];}
    void setcol(Matrix33 ma,Matrix33 mb,Matrix33 mc)
        {element[0][0]=ma[0][0];element[0][1]=mb[0][1];element[0][2]=mc[0][2];
        element[1][0]=ma[1][0];element[1][1]=mb[1][1];element[1][2]=mc[1][2];
        element[2][0]=ma[2][0];element[2][1]=mb[2][1];element[2][2]=mc[2][2];}
    void copytoarray(double e[]) const
        {e[0]=element[0][0];e[1]=element[0][1];e[2]=element[0][2];
        e[3]=element[1][0];e[4]=element[1][1];e[5]=element[1][2];
        e[6]=element[2][0];e[7]=element[2][1];e[8]=element[2][2];}
    const double * operator [] (int i) const
    { if(i>2||i<0) FATAL("wrong index Matrix["<<i<<"]");
      return element[i]; }
    double * operator [] (int i) 
    { if(i>2||i<0) FATAL("wrong index Matrix["<<i<<"]");
      return element[i]; }
    //ostream
    friend LOStream & operator <<(LOStream &os, const Matrix33 &m)
    {//formatted output
        os.Printf("%18.10e %18.10e %18.10e\n"
                "%18.10e %18.10e %18.10e\n"
                "%18.10e %18.10e %18.10e",
                m[0][0],m[0][1],m[0][2],
                m[1][0],m[1][1],m[1][2],
                m[2][0],m[2][1],m[2][2]);
        return os;
    }
    friend LIStream & operator >>(LIStream &is, Matrix33 &m)
    {
        DUMP("read matrix");
        is >> (double &)m[0][0] >> (double &)m[0][1] >> (double &)m[0][2]
           >> (double &)m[1][0] >> (double &)m[1][1] >> (double &)m[1][2]
           >> (double &)m[2][0] >> (double &)m[2][1] >> (double &)m[2][2];
        DUMP("m=[\n"<<m<<"\n");
        return is;
    }
    Matrix33 operator + (double n) const
    {
        Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]+n;
        return m;
    }
    Matrix33 operator - (double n) const
    {
        Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]-n;
        return m;
    }
    Matrix33 operator - () const
    {
        Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=-element[i][j];
        return m;
    }
    Matrix33 operator * (double n) const
    {
        Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]*n;
        return m;
    }
    Matrix33 operator / (double n) const
    {
        Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]/n;
        return m;
    }
    Matrix33 operator + (Matrix33 h) const
    {
        Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]+h[i][j];
        return m;
    }
    Matrix33 operator - (Matrix33 h) const
    {
        Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]-h[i][j];
        return m;
    }
    Matrix33 & operator += (const Matrix33 &h)
    {
        for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        element[i][j]+=h[i][j];
        return (*this);
    }
    double trace() const
    {
        return (element[0][0]+element[1][1]+element[2][2]);
    }
    void operator *= (double n)
    {
        element[0][0]*=n;
        element[0][1]*=n;
        element[0][2]*=n;
        element[1][0]*=n;
        element[1][1]*=n;
        element[1][2]*=n;
        element[2][0]*=n;
        element[2][1]*=n;
        element[2][2]*=n;
    }
    void operator /= (double n)
    {
        element[0][0]/=n;
        element[0][1]/=n;
        element[0][2]/=n;
        element[1][0]/=n;
        element[1][1]/=n;
        element[1][2]/=n;
        element[2][0]/=n;
        element[2][1]/=n;
        element[2][2]/=n;
    }
    Vector3 operator * (const Vector3 &a) const
    {
        Vector3 b;
        b.x = element[0][0] * a.x + element[0][1] * a.y + element[0][2] * a.z;
        b.y = element[1][0] * a.x + element[1][1] * a.y + element[1][2] * a.z;
        b.z = element[2][0] * a.x + element[2][1] * a.y + element[2][2] * a.z;
        return b;
    }
    void multiply (const Vector3 &a, Vector3 &b) const
    {
        b.x = element[0][0] * a.x + element[0][1] * a.y + element[0][2] * a.z;
        b.y = element[1][0] * a.x + element[1][1] * a.y + element[1][2] * a.z;
        b.z = element[2][0] * a.x + element[2][1] * a.y + element[2][2] * a.z;
    }
    
    Matrix33 operator * (const Matrix33 &m) const
    {
        Matrix33 n;
        n[0][0]=element[0][0]*m[0][0]+element[0][1]*m[1][0]+element[0][2]*m[2][0];
        n[0][1]=element[0][0]*m[0][1]+element[0][1]*m[1][1]+element[0][2]*m[2][1];
        n[0][2]=element[0][0]*m[0][2]+element[0][1]*m[1][2]+element[0][2]*m[2][2];
        n[1][0]=element[1][0]*m[0][0]+element[1][1]*m[1][0]+element[1][2]*m[2][0];
        n[1][1]=element[1][0]*m[0][1]+element[1][1]*m[1][1]+element[1][2]*m[2][1];
        n[1][2]=element[1][0]*m[0][2]+element[1][1]*m[1][2]+element[1][2]*m[2][2];
        n[2][0]=element[2][0]*m[0][0]+element[2][1]*m[1][0]+element[2][2]*m[2][0];
        n[2][1]=element[2][0]*m[0][1]+element[2][1]*m[1][1]+element[2][2]*m[2][1];
        n[2][2]=element[2][0]*m[0][2]+element[2][1]*m[1][2]+element[2][2]*m[2][2];
        return n;
    }
    Matrix33 & addnvv(double n,Vector3 &a,Vector3 &b)
    {
        element[0][0]+=n*a.x*b.x;
        element[0][1]+=n*a.x*b.y;
        element[0][2]+=n*a.x*b.z;
        element[1][0]+=n*a.y*b.x;
        element[1][1]+=n*a.y*b.y;
        element[1][2]+=n*a.y*b.z;
        element[2][0]+=n*a.z*b.x;
        element[2][1]+=n*a.z*b.y;
        element[2][2]+=n*a.z*b.z;
        return (*this);
    }
    Matrix33 & dyad(Vector3 &a,Vector3 &b)
    {
        element[0][0]=a.x*b.x;
        element[0][1]=a.x*b.y;
        element[0][2]=a.x*b.z;
        element[1][0]=a.y*b.x;
        element[1][1]=a.y*b.y;
        element[1][2]=a.y*b.z;
        element[2][0]=a.z*b.x;
        element[2][1]=a.z*b.y;
        element[2][2]=a.z*b.z;
        return (*this);
    }
    Matrix33 adj() const
    {
        Matrix33 hadj;
        hadj[0][0]=element[1][1]*element[2][2]-element[1][2]*element[2][1];
        hadj[1][1]=element[2][2]*element[0][0]-element[2][0]*element[0][2];
        hadj[2][2]=element[0][0]*element[1][1]-element[0][1]*element[1][0];
        hadj[0][1]=element[1][2]*element[2][0]-element[1][0]*element[2][2];
        hadj[1][2]=element[2][0]*element[0][1]-element[2][1]*element[0][0];
        hadj[2][0]=element[0][1]*element[1][2]-element[0][2]*element[1][1];
        hadj[0][2]=element[1][0]*element[2][1]-element[2][0]*element[1][1];
        hadj[1][0]=element[2][1]*element[0][2]-element[0][1]*element[2][2];
        hadj[2][1]=element[0][2]*element[1][0]-element[1][2]*element[0][0];
        return hadj;
    }
    Matrix33 tran() const
    {
        Matrix33 htrn;
        htrn[0][0]=element[0][0];htrn[1][1]=element[1][1];htrn[2][2]=element[2][2];
        htrn[0][1]=element[1][0];htrn[1][2]=element[2][1];htrn[2][0]=element[0][2];
        htrn[1][0]=element[0][1];htrn[2][1]=element[1][2];htrn[0][2]=element[2][0];
        return htrn;
    }
    Matrix33 inv() const
    {
        Matrix33 hinv;double c;
        hinv=adj();
        c=element[0][0]*hinv[0][0]
            +element[0][1]*hinv[0][1]+element[0][2]*hinv[0][2];
        hinv=hinv.tran();
        hinv[0][0]/=c;hinv[0][1]/=c;hinv[0][2]/=c;
        hinv[1][0]/=c;hinv[1][1]/=c;hinv[1][2]/=c;
        hinv[2][0]/=c;hinv[2][1]/=c;hinv[2][2]/=c;
        return hinv;
    }
    Matrix33 reorient() const
    {
        double a, b, c, alpha, beta, gamma;
        double b1,b2,c1,c2,c3;
        Vector3 va, vb, vc;   Matrix33 h0;

        va.set(element[0][0],element[1][0],element[2][0]);
        vb.set(element[0][1],element[1][1],element[2][1]);
        vc.set(element[0][2],element[1][2],element[2][2]);
        a=va.norm();
        b=vb.norm();
        c=vc.norm();
        alpha=acos(vb*vc/b/c);
        beta=acos(vc*va/c/a);
        gamma=acos(va*vb/a/b);
        b1=b*cos(gamma);
        b2=b*sin(gamma);
        c1=c*cos(beta);
        c2=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
        c3=sqrt(c*c-c1*c1-c2*c2);
        h0.set(a,b1,c1,
               0,b2,c2,
               0,0,c3);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                if(fabs(h0[i][j])<1e-11) h0[i][j]=0.0;
        
        return h0;
        
    }
    double det() const
    {
        Matrix33 hadj;double c;
        hadj=adj();
        c=element[0][0]*hadj[0][0]
            +element[0][1]*hadj[0][1]+element[0][2]*hadj[0][2];
        return c;
    }
    Vector3 height() const
    {
        Vector3 h; Vector3 r;

        /*
        double V; Matrix33 hadj;           
        V=det();
        hadj=adj();
        r.set(hadj[0][0],hadj[1][0],hadj[2][0]);
        h[0]=V/r.norm();
        r.set(hadj[0][1],hadj[1][1],hadj[2][1]);
        h[1]=V/r.norm();
        r.set(hadj[0][2],hadj[1][2],hadj[2][2]);
        h[2]=V/r.norm();
        */
        
        Matrix33 hinv;
        hinv=inv();
        /* INFO("hinv="<<hinv); */
        r.set(hinv[0][0],hinv[0][1],hinv[0][2]); h[0]=1/r.norm();
        r.set(hinv[1][0],hinv[1][1],hinv[1][2]); h[1]=1/r.norm();
        r.set(hinv[2][0],hinv[2][1],hinv[2][2]); h[2]=1/r.norm();
        

        return h;
    }
    void clear()
    {
        element[0][0]=element[0][1]=element[0][2]=
            element[1][0]=element[1][1]=element[1][2]=
            element[2][0]=element[2][1]=element[2][2]=0;
    }
    void eye()
    {
        element[0][1]=element[0][2]=
            element[1][0]=element[1][2]=
            element[2][0]=element[2][1]=0;
        element[0][0]=element[1][1]=element[2][2]=1;
    }
};


class UnitCell
{
/* cannot change this static memory into dynamic allocation */
#define MAXBASISNUM 20000
//#define MAXBASISNUM 1000
public:
    Matrix33 h;
    int n;
    Vector3 basis[MAXBASISNUM];
    int species[MAXBASISNUM];
    
    UnitCell():n(1)
    {
        h.eye();
        basis[0].x=basis[0].y=basis[0].z=0;
        species[0]=0;
    }
    UnitCell(int m,const double *b)
    {
        h.eye();
        n=m;
        if(n>=MAXBASISNUM)
            FATAL("UnitCell n="<<n<<" exceed limit"<<MAXBASISNUM);
        for(int i=0;i<n;i++)
        {
            basis[i].clear();
            basis[i].x=b[3*i];
            basis[i].y=b[3*i+1];
            basis[i].z=b[3*i+2];
            species[i]=0;
        }
    }
    UnitCell(int m,const double *b,const int *s)
    {
        h.eye();
        n=m;
        if(n>=MAXBASISNUM)
            FATAL("UnitCell n="<<n<<" exceed limit"<<MAXBASISNUM);
        for(int i=0;i<n;i++)
        {
            basis[i].clear();
            basis[i].x=b[3*i];
            basis[i].y=b[3*i+1];
            basis[i].z=b[3*i+2];
            species[i]=s[i];
        }
    }
    void set(int m,const double *b)
    {
        h.eye();
        n=m;
        if(n>=MAXBASISNUM)
            FATAL("UnitCell n="<<n<<" exceed limit"<<MAXBASISNUM);
        for(int i=0;i<n;i++)
        {
            basis[i].clear();
            basis[i].x=b[3*i];
            basis[i].y=b[3*i+1];
            basis[i].z=b[3*i+2];
            species[i]=0;
        }
    }
    void set(int m,const double *b,const int *s)
    {
        h.eye();
        n=m;
        if(n>=MAXBASISNUM)
            FATAL("UnitCell n="<<n<<" exceed limit"<<MAXBASISNUM);
        for(int i=0;i<n;i++)
        {
            basis[i].clear();
            basis[i].x=b[3*i];
            basis[i].y=b[3*i+1];
            basis[i].z=b[3*i+2];
            species[i]=s[i];
        }
    }
    UnitCell operator * (const Matrix33 &m) const
    {
        Vector3 r[8];//convexity of unit cell, 8 extreme points
        Vector3 min,max;
        Vector3 s,offset,ns;
        Matrix33 hinv;
        int i,j,k,p,q;
        
        UnitCell uc;
        uc=(*this);
        uc.h=h*m;
        hinv=m.inv();

        r[0].set(0,0,0);
        r[1].set(m[0][0],m[1][0],m[2][0]);
        r[2].set(m[0][1],m[1][1],m[2][1]);
        r[3].set(m[0][2],m[1][2],m[2][2]);
        r[4]=r[1]+r[2];
        r[5]=r[2]+r[3];
        r[6]=r[3]+r[1];
        r[7]=r[3]+r[1]+r[2];

        min=r[0];max=r[0];
        for(i=1;i<8;i++)
        {
            if(min.x>r[i].x) min.x=r[i].x;
            if(min.y>r[i].y) min.y=r[i].y;
            if(min.z>r[i].z) min.z=r[i].z;
            if(max.x<r[i].x) max.x=r[i].x;
            if(max.y<r[i].y) max.y=r[i].y;
            if(max.z<r[i].z) max.z=r[i].z;
        }

        DUMP("min="<<min<<"\nmax="<<max);
        q=0;
        for(i=(int)floor(min.x);i<=(int)ceil(max.x);i++)
            for(j=(int)floor(min.y);j<=(int)ceil(max.y);j++)
                for(k=(int)floor(min.z);k<=(int)ceil(max.z);k++)
                {
                    if(n>=MAXBASISNUM)
                        FATAL("UnitCell n="<<n<<" exceed limit"<<MAXBASISNUM);
                    for(p=0;p<n;p++)
                    {
                        offset.set(i,j,k);
                        s=basis[p]+offset;
                        offset.set(0.5,0.5,0.5);
                        ns=hinv*s;
#define margin 1e-8
#define within01(x) (((x)>-margin)&&((x)<1-margin))
                        if((within01(ns.x))&&(within01(ns.y))&&(within01(ns.z)))
                        {
                            uc.basis[q]=ns;
                            uc.species[q]=species[p];
                            if(q>=MAXBASISNUM)
                                FATAL("UnitCell q="<<q<<" exceed limit"<<MAXBASISNUM);
                            q++;
                            DUMP("UnitCell: q="<<q<<"\ns="<<s<<"\nns="<<ns);
                        }
                    }
                }
        uc.n=q;
        return uc;
    }
    friend LOStream & operator <<(LOStream &os, const UnitCell &uc)
    {//formatted output
        os.Printf("h=[\n%10g %10g %10g\n"
                "%10g %10g %10g\n"
                "%10g %10g %10g];\n",
                uc.h[0][0],uc.h[0][1],uc.h[0][2],
                uc.h[1][0],uc.h[1][1],uc.h[1][2],
                uc.h[2][0],uc.h[2][1],uc.h[2][2]);
        os<<"basis=[\n";
        for(int i=0;i<uc.n;i++)
            os.Printf("%10g %10g %10g %2d   %6d\n",
                      uc.basis[i].x,uc.basis[i].y,uc.basis[i].z,uc.species[i],i);
        os<<"];";
        return os;
    }
};

#endif // _LINALG3_H

