/*
  linalg3.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sun Dec  3 16:06:53 2006

  FUNCTION  :  Linear Algebra package for 3 dimension array and matrix

  Featuring :  1. 3x3 Matrix
               2. 3 Vector
               3. UnitCell
*/

#ifndef _LINALG3_CUDA_H
#define _LINALG3_CUDA_H

#include <math.h>

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

class G_Vector3
{
public:
 double x, y, z;
 CUDA_CALLABLE_MEMBER   G_Vector3():x(0),y(0),z(0){};
 CUDA_CALLABLE_MEMBER   G_Vector3(const Vector3& a){x=a.x; y=a.y;z=a.z;};
 CUDA_CALLABLE_MEMBER   double operator [] (int i) const { return (&x)[i]; }
 CUDA_CALLABLE_MEMBER   double & operator [] (int i) { return (&x)[i]; }
 CUDA_CALLABLE_MEMBER   G_Vector3(double a,double b,double c):x(a),y(b),z(c){};
 CUDA_CALLABLE_MEMBER   G_Vector3(const double a[]){x=a[0];y=a[1];z=a[2];}
 CUDA_CALLABLE_MEMBER   
 CUDA_CALLABLE_MEMBER   G_Vector3 operator + (const G_Vector3 &a) const
 CUDA_CALLABLE_MEMBER   { G_Vector3 c; c.x=x+a.x;c.y=y+a.y;c.z=z+a.z; return c; }
 CUDA_CALLABLE_MEMBER   G_Vector3 operator - (G_Vector3 &a) const
 CUDA_CALLABLE_MEMBER   { G_Vector3 c; c.x=x-a.x;c.y=y-a.y;c.z=z-a.z; return c; }
 CUDA_CALLABLE_MEMBER   G_Vector3 operator * (double b) const
 CUDA_CALLABLE_MEMBER   { G_Vector3 c; c.x=x*b;c.y=y*b;c.z=z*b; return c; }
 CUDA_CALLABLE_MEMBER   G_Vector3 operator / (double b) const
 CUDA_CALLABLE_MEMBER   { G_Vector3 c; c.x=x/b;c.y=y/b;c.z=z/b; return c; }

 CUDA_CALLABLE_MEMBER   void add     (G_Vector3 &a, G_Vector3& b){x=a.x+b.x; y=a.y+b.y; z=a.z+b.z;}
 CUDA_CALLABLE_MEMBER   void subtract(G_Vector3 &a, G_Vector3& b){x=a.x-b.x; y=a.y-b.y; z=a.z-b.z;}
 CUDA_CALLABLE_MEMBER   void addnv(double a,const G_Vector3& b){x+=a*b.x;y+=a*b.y;z+=a*b.z;}
 CUDA_CALLABLE_MEMBER   void invert(){x=1./x;y=1./y;z=1./z;}
 CUDA_CALLABLE_MEMBER   void subint() {x=x-rint(x);y=y-rint(y);z=z-rint(z);}
 CUDA_CALLABLE_MEMBER   double norm2() const {return x*x+y*y+z*z; }
 CUDA_CALLABLE_MEMBER   double norm() const  {return sqrt(x*x+y*y+z*z); }
 CUDA_CALLABLE_MEMBER   void clear() {x=y=z=0;}
 CUDA_CALLABLE_MEMBER   double operator * (const G_Vector3 &a) const
                        { return (x*a.x+y*a.y+z*a.z); }
    
#if 0 /* return value */
    G_Vector3 & operator += (double b) { x+=b;y+=b;z+=b; return (*this); }
    G_Vector3 & operator -= (double b) { x-=b;y-=b;z-=b; return (*this); }
    G_Vector3 & operator *= (double b) { x*=b;y*=b;z*=b; return (*this); }
    G_Vector3 & operator /= (double b) { x/=b;y/=b;z/=b; return (*this); }
    G_Vector3 & operator += (const G_Vector3 &a) { x+=a.x;y+=a.y;z+=a.z; return (*this); }
    G_Vector3 & operator -= (const G_Vector3 &a) { x-=a.x;y-=a.y;z-=a.z; return (*this); }
    G_Vector3 & operator *= (const G_Vector3 &a) { x*=a.x;y*=a.y;z*=a.z; return (*this); }
    G_Vector3 & operator /= (const G_Vector3 &a) { x/=a.x;y/=a.y;z/=a.z; return (*this); }
#else /* do not return value */
 CUDA_CALLABLE_MEMBER   void operator += (double b) { x+=b;y+=b;z+=b; }
 CUDA_CALLABLE_MEMBER   void operator -= (double b) { x-=b;y-=b;z-=b; }
 CUDA_CALLABLE_MEMBER   void operator *= (double b) { x*=b;y*=b;z*=b; }
 CUDA_CALLABLE_MEMBER   void operator /= (double b) { x/=b;y/=b;z/=b; }
 CUDA_CALLABLE_MEMBER   void operator += (const G_Vector3 &a) { x+=a.x;y+=a.y;z+=a.z; }
 CUDA_CALLABLE_MEMBER   void operator -= (const G_Vector3 &a) { x-=a.x;y-=a.y;z-=a.z; }
 CUDA_CALLABLE_MEMBER   void operator *= (const G_Vector3 &a) { x*=a.x;y*=a.y;z*=a.z; }
 CUDA_CALLABLE_MEMBER   void operator /= (const G_Vector3 &a) { x/=a.x;y/=a.y;z/=a.z; }
#endif

 CUDA_CALLABLE_MEMBER   bool operator == (const G_Vector3 &b) {
   if(fabs(x- b.x)<1e-9&&fabs(y- b.y)<1e-9&&fabs(z- b.z)<1e-9)
     return true; 
   else 
     return false; 
 }
    
 CUDA_CALLABLE_MEMBER   void set(const double a,const double b,const double c){x=a;y=b;z=c;}
 CUDA_CALLABLE_MEMBER   void set(const double a[]){x=a[0];y=a[1];z=a[2];}
 CUDA_CALLABLE_MEMBER   void add(const double a,const double b,const double c){x+=a;y+=b;z+=c;}
 CUDA_CALLABLE_MEMBER   void copytoarray(double a[]) const {a[0]=x;a[1]=y;a[2]=z;}
 CUDA_CALLABLE_MEMBER   G_Vector3 sq() const { G_Vector3 c; c.x=x*x;c.y=y*y;c.z=z*z; return c;}
 CUDA_CALLABLE_MEMBER   G_Vector3 sqroot() const { G_Vector3 c; c.x=(double)sqrt(x);c.y=sqrt(y);c.z=sqrt(z); return c;}
    
  CUDA_CALLABLE_MEMBER  G_Vector3 & orth(G_Vector3 &a)
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
 CUDA_CALLABLE_MEMBER   G_Vector3 & proj(G_Vector3 &a)
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
 CUDA_CALLABLE_MEMBER   friend double dot(G_Vector3 &a, G_Vector3 &b)
    {
        double c;
        c=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        return c;
    }
 CUDA_CALLABLE_MEMBER   friend G_Vector3 cross(G_Vector3 &a, G_Vector3 &b)
    {
        G_Vector3 c;
        c[0]=a[1]*b[2]-a[2]*b[1];
        c[1]=a[2]*b[0]-a[0]*b[2];
        c[2]=a[0]*b[1]-a[1]*b[0];
        return c;
    }
 CUDA_CALLABLE_MEMBER   friend G_Vector3 schmidt(G_Vector3 &a, G_Vector3 &b)
    {
#ifndef TINYNORM2
#define TINYNORM2 1e-10
#endif
        G_Vector3 c,a1;
        a1=a;
        c=a1;c.orth(b);
        if(c.norm2()<TINYNORM2) a1[0]+=1;
        c=a1;c.orth(b);
        if(c.norm2()<TINYNORM2) a1[1]+=1;
        c=a1;c.orth(b);
        if(c.norm2()<TINYNORM2)
        {
           // FATAL("Schmidt("<<a<<", "<<b<<"): "
            //      <<"too small! c="<<c);
        }
#undef TINYNORM
        return c/c.norm();
    }
 CUDA_CALLABLE_MEMBER   friend LOStream & operator <<(LOStream &os, const G_Vector3 &a)
    {//formatted output
        //os.Printf("%18.10e %18.10e %18.10e",a.x,a.y,a.z);
        return os;
    }
 CUDA_CALLABLE_MEMBER   friend LIStream & operator >>(LIStream &is, G_Vector3 &a)
    {
        is >> a.x >> a.y >> a.z;
        return is;
    }

};

class G_Matrix33
{
public:
 double element[3][3];
public:
 CUDA_CALLABLE_MEMBER   G_Matrix33(){clear();}
 CUDA_CALLABLE_MEMBER   G_Matrix33(const Matrix33 e)
        {element[0][0]=e.element[0][0];element[0][1]=e.element[0][1];element[0][2]=e.element[0][2];
         element[1][0]=e.element[1][0];element[1][1]=e.element[1][1];element[1][2]=e.element[1][2];
         element[2][0]=e.element[2][0];element[2][1]=e.element[2][1];element[2][2]=e.element[2][2];}
 CUDA_CALLABLE_MEMBER   G_Matrix33(double e1,double e2,double e3,
             double e4,double e5,double e6,
             double e7,double e8,double e9)
        {element[0][0]=e1;element[0][1]=e2;element[0][2]=e3;
        element[1][0]=e4;element[1][1]=e5;element[1][2]=e6;
        element[2][0]=e7;element[2][1]=e8;element[2][2]=e9;}
 CUDA_CALLABLE_MEMBER   G_Matrix33(double e[9])
        {element[0][0]=e[0];element[0][1]=e[1];element[0][2]=e[2];
        element[1][0]=e[3];element[1][1]=e[4];element[1][2]=e[5];
        element[2][0]=e[6];element[2][1]=e[7];element[2][2]=e[8];}
 CUDA_CALLABLE_MEMBER   void set(const double e1,const double e2,const double e3,
             const double e4,const double e5,const double e6,
             const double e7,const double e8,const double e9)
        {element[0][0]=e1;element[0][1]=e2;element[0][2]=e3;
        element[1][0]=e4;element[1][1]=e5;element[1][2]=e6;
        element[2][0]=e7;element[2][1]=e8;element[2][2]=e9;}
 CUDA_CALLABLE_MEMBER   void set(const double e[])
        {element[0][0]=e[0];element[0][1]=e[1];element[0][2]=e[2];
        element[1][0]=e[3];element[1][1]=e[4];element[1][2]=e[5];
        element[2][0]=e[6];element[2][1]=e[7];element[2][2]=e[8];}
 CUDA_CALLABLE_MEMBER   void setcol(G_Vector3 &a,G_Vector3 &b,G_Vector3 &c)
        {element[0][0]=a[0];element[0][1]=b[0];element[0][2]=c[0];
        element[1][0]=a[1];element[1][1]=b[1];element[1][2]=c[1];
        element[2][0]=a[2];element[2][1]=b[2];element[2][2]=c[2];}
 CUDA_CALLABLE_MEMBER   void setcol(double a[],double b[],double c[])
        {element[0][0]=a[0];element[0][1]=b[0];element[0][2]=c[0];
        element[1][0]=a[1];element[1][1]=b[1];element[1][2]=c[1];
        element[2][0]=a[2];element[2][1]=b[2];element[2][2]=c[2];}
 CUDA_CALLABLE_MEMBER   void setcol(G_Matrix33 ma,G_Matrix33 mb,G_Matrix33 mc)
        {element[0][0]=ma[0][0];element[0][1]=mb[0][1];element[0][2]=mc[0][2];
        element[1][0]=ma[1][0];element[1][1]=mb[1][1];element[1][2]=mc[1][2];
        element[2][0]=ma[2][0];element[2][1]=mb[2][1];element[2][2]=mc[2][2];}
 CUDA_CALLABLE_MEMBER   void copytoarray(double e[]) const
        {e[0]=element[0][0];e[1]=element[0][1];e[2]=element[0][2];
        e[3]=element[1][0];e[4]=element[1][1];e[5]=element[1][2];
        e[6]=element[2][0];e[7]=element[2][1];e[8]=element[2][2];}
 CUDA_CALLABLE_MEMBER   const double * operator [] (int i) const
    { //if(i>2||i<0) FATAL("wrong index Matrix["<<i<<"]");
      return element[i]; }
 CUDA_CALLABLE_MEMBER   double * operator [] (int i) 
    { //if(i>2||i<0) FATAL("wrong index Matrix["<<i<<"]");
      return element[i]; }
    //ostream
 CUDA_CALLABLE_MEMBER   friend LOStream & operator <<(LOStream &os, const G_Matrix33 &m)
    {//formatted output
      //  os.Printf("%18.10e %18.10e %18.10e\n"
        //        "%18.10e %18.10e %18.10e\n"
        //        "%18.10e %18.10e %18.10e",
        //        m[0][0],m[0][1],m[0][2],
        //        m[1][0],m[1][1],m[1][2],
        //        m[2][0],m[2][1],m[2][2]);
        return os;
    }
 CUDA_CALLABLE_MEMBER   friend LIStream & operator >>(LIStream &is, G_Matrix33 &m)
    {
        DUMP("read matrix");
        is >> (double &)m[0][0] >> (double &)m[0][1] >> (double &)m[0][2]
           >> (double &)m[1][0] >> (double &)m[1][1] >> (double &)m[1][2]
           >> (double &)m[2][0] >> (double &)m[2][1] >> (double &)m[2][2];
        DUMP("m=[\n"<<m<<"\n");
        return is;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 operator + (double n) const
    {
        G_Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]+n;
        return m;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 operator - (double n) const
    {
        G_Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]-n;
        return m;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 operator - () const
    {
        G_Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=-element[i][j];
        return m;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 operator * (double n) const
    {
        G_Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]*n;
        return m;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 operator / (double n) const
    {
        G_Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]/n;
        return m;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 operator + (G_Matrix33 h) const
    {
        G_Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]+h[i][j];
        return m;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 operator - (G_Matrix33 h) const
    {
        G_Matrix33 m;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                m[i][j]=element[i][j]-h[i][j];
        return m;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 & operator += (const G_Matrix33 &h)
    {
        for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        element[i][j]+=h[i][j];
        return (*this);
    }
 CUDA_CALLABLE_MEMBER   double trace() const
    {
        return (element[0][0]+element[1][1]+element[2][2]);
    }
 CUDA_CALLABLE_MEMBER   void operator *= (double n)
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
 CUDA_CALLABLE_MEMBER   void operator /= (double n)
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
 CUDA_CALLABLE_MEMBER   G_Vector3 operator * (const G_Vector3 &a) const
    {
        G_Vector3 b;
        b.x = element[0][0] * a.x + element[0][1] * a.y + element[0][2] * a.z;
        b.y = element[1][0] * a.x + element[1][1] * a.y + element[1][2] * a.z;
        b.z = element[2][0] * a.x + element[2][1] * a.y + element[2][2] * a.z;
        return b;
    }
 CUDA_CALLABLE_MEMBER   void multiply (const G_Vector3 &a, G_Vector3 &b) const
    {
        b.x = element[0][0] * a.x + element[0][1] * a.y + element[0][2] * a.z;
        b.y = element[1][0] * a.x + element[1][1] * a.y + element[1][2] * a.z;
        b.z = element[2][0] * a.x + element[2][1] * a.y + element[2][2] * a.z;
    }
    
 CUDA_CALLABLE_MEMBER   G_Matrix33 operator * (const G_Matrix33 &m) const
    {
        G_Matrix33 n;
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
 CUDA_CALLABLE_MEMBER   bool operator == (const G_Matrix33 &b) {
	for(int i = 0;i<3;i++) for(int j = 0;j<3;j++) {
          if(fabs(element[i][j] -b.element[i][j])>1e-15) 
            return false; 
	}
        return true; 
 }

 CUDA_CALLABLE_MEMBER   G_Matrix33 & addnvv(double n,G_Vector3 &a,G_Vector3 &b)
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
 CUDA_CALLABLE_MEMBER   G_Matrix33 & dyad(G_Vector3 &a,G_Vector3 &b)
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
 CUDA_CALLABLE_MEMBER   G_Matrix33 adj() const
    {
        G_Matrix33 hadj;
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
 CUDA_CALLABLE_MEMBER   G_Matrix33 tran() const
    {
        G_Matrix33 htrn;
        htrn[0][0]=element[0][0];htrn[1][1]=element[1][1];htrn[2][2]=element[2][2];
        htrn[0][1]=element[1][0];htrn[1][2]=element[2][1];htrn[2][0]=element[0][2];
        htrn[1][0]=element[0][1];htrn[2][1]=element[1][2];htrn[0][2]=element[2][0];
        return htrn;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 inv() const
    {
        G_Matrix33 hinv;double c;
        hinv=adj();
        c=element[0][0]*hinv[0][0]
            +element[0][1]*hinv[0][1]+element[0][2]*hinv[0][2];
        hinv=hinv.tran();
        hinv[0][0]/=c;hinv[0][1]/=c;hinv[0][2]/=c;
        hinv[1][0]/=c;hinv[1][1]/=c;hinv[1][2]/=c;
        hinv[2][0]/=c;hinv[2][1]/=c;hinv[2][2]/=c;
        return hinv;
    }
 CUDA_CALLABLE_MEMBER   G_Matrix33 reorient() const
    {
        double a, b, c, alpha, beta, gamma;
        double b1,b2,c1,c2,c3;
        G_Vector3 va, vb, vc;   G_Matrix33 h0;

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
 CUDA_CALLABLE_MEMBER   double det() const
    {
        G_Matrix33 hadj;double c;
        hadj=adj();
        c=element[0][0]*hadj[0][0]
            +element[0][1]*hadj[0][1]+element[0][2]*hadj[0][2];
        return c;
    }
 CUDA_CALLABLE_MEMBER   G_Vector3 height() const
    {
        G_Vector3 h; G_Vector3 r;

        /*
        double V; G_Matrix33 hadj;           
        V=det();
        hadj=adj();
        r.set(hadj[0][0],hadj[1][0],hadj[2][0]);
        h[0]=V/r.norm();
        r.set(hadj[0][1],hadj[1][1],hadj[2][1]);
        h[1]=V/r.norm();
        r.set(hadj[0][2],hadj[1][2],hadj[2][2]);
        h[2]=V/r.norm();
        */
        
        G_Matrix33 hinv;
        hinv=inv();
        /* INFO("hinv="<<hinv); */
        r.set(hinv[0][0],hinv[0][1],hinv[0][2]); h[0]=1/r.norm();
        r.set(hinv[1][0],hinv[1][1],hinv[1][2]); h[1]=1/r.norm();
        r.set(hinv[2][0],hinv[2][1],hinv[2][2]); h[2]=1/r.norm();
        

        return h;
    }
 CUDA_CALLABLE_MEMBER   void clear()
    {
        element[0][0]=element[0][1]=element[0][2]=
            element[1][0]=element[1][1]=element[1][2]=
            element[2][0]=element[2][1]=element[2][2]=0;
    }
 CUDA_CALLABLE_MEMBER   void eye()
    {
        element[0][1]=element[0][2]=
            element[1][0]=element[1][2]=
            element[2][0]=element[2][1]=0;
        element[0][0]=element[1][1]=element[2][2]=1;
    }
};



#endif // _LINALG3_CUDA_H

