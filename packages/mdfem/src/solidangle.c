/* this file contains subroutines to compute solid angles of triangles, polygons, and ellipse */
#include <math.h>
#include <stdio.h>

/*-------------------------------------------------------------------------
 *
 *      Function:       SolidAngleTriangle
 *      Description:    Calculate the solid angle of a triangle.
 *
 *      Arguments:
 *              x1,y1,z1     First  vertex of the triangle
 *              x2,y2,z2     Second vertex of the triangle
 *              x3,y3,z3     Third  vertex of the triangle
 *              xp,yp,zp     Point at which solid angle is computed
 *      Return:
 *              Omega        solid angle
 *-----------------------------------------------------------------------*/
double SolidAngleTriangle(double x1, double y1, double z1,
                        double x2, double y2, double z2,
                        double x3, double y3, double z3, 
                        double xp, double yp, double zp)
{
    double xs1,ys1,zs1,xs2,ys2,zs2,xs3,ys3,zs3;
    double RR1,RR2,RR3,dotR1R2,dotR2R3,dotR3R1;
    double numer,denom;
    double Omega;
    
    xs1 = x1 - xp; ys1 = y1 - yp; zs1 = z1 - zp;
    xs2 = x2 - xp; ys2 = y2 - yp; zs2 = z2 - zp;
    xs3 = x3 - xp; ys3 = y3 - yp; zs3 = z3 - zp;

    RR1 = sqrt(xs1*xs1 + ys1*ys1 + zs1*zs1);
    RR2 = sqrt(xs2*xs2 + ys2*ys2 + zs2*zs2);
    RR3 = sqrt(xs3*xs3 + ys3*ys3 + zs3*zs3);
    
    numer = xs1*ys2*zs3 + ys1*zs2*xs3 + zs1*xs2*ys3
           -xs1*zs2*ys3 - ys1*xs2*zs3 - zs1*ys2*xs3;

    dotR1R2 = xs1*xs2 + ys1*ys2 + zs1*zs2;
    dotR2R3 = xs2*xs3 + ys2*ys3 + zs2*zs3;
    dotR3R1 = xs3*xs1 + ys3*ys1 + zs3*zs1;
    
    denom = RR1*RR2*RR3 + dotR1R2*RR3 + dotR2R3*RR1 + dotR3R1*RR2;
  
    Omega = -2 * atan2(numer, denom);

    return Omega;
}

Vector3 get_f_Triangle(double x1, double y1, double z1,
                         double x2, double y2, double z2,
                         double x3, double y3, double z3, 
                         double xp, double yp, double zp, Vector3 *burg)
{
    double xs1,ys1,zs1,xs2,ys2,zs2,xs3,ys3,zs3;
    double norm_A, norm_B, coeff_f;
    class Vector3 R1, R2, R3, df, f, t_AB, lambda_A, lambda_B;
    
    xs1 = x1 - xp; ys1 = y1 - yp; zs1 = z1 - zp;
    xs2 = x2 - xp; ys2 = y2 - yp; zs2 = z2 - zp;
    xs3 = x3 - xp; ys3 = y3 - yp; zs3 = z3 - zp;

    R1.set(xs1, ys1, zs1);
    R2.set(xs2, ys2, zs2);
    R3.set(xs3, ys3, zs3);

    f.set(0,0,0);
    
#define get_f(A,B)  protect(\
        /* unit tangent */                                       \
        t_AB = B - A;                                            \
        t_AB /= sqrt(t_AB.x*t_AB.x+t_AB.y*t_AB.y+t_AB.z*t_AB.z); \
        /* unit vectors along R1 and R2 */      \
        norm_A = sqrt(A.x*A.x+A.y*A.y+A.z*A.z); lambda_A = A/norm_A; \
        norm_B = sqrt(B.x*B.x+B.y*B.y+B.z*B.z); lambda_B = B/norm_B; \
        coeff_f=log(norm_B/norm_A*(1+dot(lambda_B,t_AB))/(1+dot(lambda_A,t_AB))); \
        df=cross(*burg,t_AB)*coeff_f;             \
        )
  
    get_f(R1,R2); f = df;
    get_f(R2,R3); f = f + df;
    get_f(R3,R1); f = f + df;
    
//    INFO_Printf("In get_f_Triangle, b=[%f %f %f]\n",(*burg).x,(*burg).y,(*burg).z);
//    INFO_Printf("f = [%f %f %f]\n",f.x,f.y,f.z);

    return f;
}

Vector3 get_g_Triangle(double x1, double y1, double z1,
                       double x2, double y2, double z2,
                       double x3, double y3, double z3, 
                       double xp, double yp, double zp, Vector3 *burg)
{
    double xs1,ys1,zs1,xs2,ys2,zs2,xs3,ys3,zs3;
    double norm_A, norm_B, coeff_g;
    class Vector3 R1, R2, R3, dg, g, lambda_A, lambda_B, cross_AB;
    
    xs1 = x1 - xp; ys1 = y1 - yp; zs1 = z1 - zp;
    xs2 = x2 - xp; ys2 = y2 - yp; zs2 = z2 - zp;
    xs3 = x3 - xp; ys3 = y3 - yp; zs3 = z3 - zp;

    R1.set(xs1, ys1, zs1);
    R2.set(xs2, ys2, zs2);
    R3.set(xs3, ys3, zs3);

    g.set(0,0,0);
    
#define get_g(A,B)  protect(\
        /* unit vectors along R1 and R2 */      \
        norm_A = sqrt(A.x*A.x+A.y*A.y+A.z*A.z); lambda_A = A/norm_A; \
        norm_B = sqrt(B.x*B.x+B.y*B.y+B.z*B.z); lambda_B = B/norm_B; \
        cross_AB = cross(lambda_A,lambda_B); \
        coeff_g=dot(*burg,cross_AB)/(1+dot(lambda_A,lambda_B)); \
        dg=(lambda_A + lambda_B)*coeff_g;        \
        )
  
    get_g(R1,R2); g = dg;
    get_g(R2,R3); g = g + dg;
    get_g(R3,R1); g = g + dg;
    
//    INFO_Printf("In get_g_Triangle, b=[%f %f %f]\n",(*burg).x,(*burg).y,(*burg).z);
//    INFO_Printf("g = [%f %f %f]\n",g.x,g.y,g.z);

    return g;
}

/*-------------------------------------------------------------------------
 *
 *      Function:       SolidAnglePolygon
 *      Description:    Calculate the solid angle of a polygon
 *
 *      Arguments:
 *              N            Number of sides of the polygon
 *              rn[3*N]      coordinates of the polygon vertices (x1,y1,z1,x2,y2,z2,...xn,yn,zn)
 *              xp,yp,zp     Point at which solid angle is computed
 *      Return:
 *              Omega        solid angle
 *-----------------------------------------------------------------------*/
double SolidAnglePolygon (int N,
                          double rn[],
                          double xp,  double yp,  double zp )
{
    int i;
    double Omega;

    if (N<3) return 0;
    
    Omega = SolidAngleTriangle(rn[0],rn[1],rn[2],
                               rn[3],rn[4],rn[5],
                               rn[6],rn[7],rn[8],
                               xp,   yp,   zp  );

    for(i=3;i<N;i++)
        Omega += SolidAngleTriangle(rn[0],rn[1],rn[2],
                                    rn[3*(i-1)],rn[3*(i-1)+1],rn[3*(i-1)+2],
                                    rn[3*i],rn[3*i+1],rn[3*i+2],
                                    xp,   yp,   zp  );

    return Omega;
}

Vector3 get_f_Polygon (int N,
                       double rn[],
                       double xp,  double yp,  double zp, Vector3 *b)
{
    int i;
    class Vector3 f;

    if (N<3) return 0;
    
    f = get_f_Triangle(rn[0],rn[1],rn[2],
                       rn[3],rn[4],rn[5],
                       rn[6],rn[7],rn[8],
                       xp,   yp,   zp,   b);

    for(i=3;i<N;i++)
        f += get_f_Triangle(rn[0],rn[1],rn[2],
                            rn[3*(i-1)],rn[3*(i-1)+1],rn[3*(i-1)+2],
                            rn[3*i],rn[3*i+1],rn[3*i+2],
                            xp,   yp,   zp,   b);
    
//    INFO_Printf("In get_f_Polygon, b=[%f %f %f]\n",(*b).x,(*b).y,(*b).z);
//    INFO_Printf("f = [%f %f %f]\n",f.x,f.y,f.z);
    
    return f;
}

Vector3 get_g_Polygon (int N,
                       double rn[],
                       double xp,  double yp,  double zp, Vector3 *b )
{
    int i;
    class Vector3 g;

    if (N<3) return 0;
    
    g = get_g_Triangle(rn[0],rn[1],rn[2],
                       rn[3],rn[4],rn[5],
                       rn[6],rn[7],rn[8],
                       xp,   yp,   zp,   b);

    for(i=3;i<N;i++)
        g += get_g_Triangle(rn[0],rn[1],rn[2],
                            rn[3*(i-1)],rn[3*(i-1)+1],rn[3*(i-1)+2],
                            rn[3*i],rn[3*i+1],rn[3*i+2],
                            xp,   yp,   zp,   b);

//    INFO_Printf("In get_g_Polygon, b=[%f %f %f]\n",(*b).x,(*b).y,(*b).z);
//    INFO_Printf("g = [%f %f %f]\n",g.x,g.y,g.z);
    
    return g;
}

/*-------------------------------------------------------------------------
 *
 *      Function:       SolidAngleEllipse
 *      Description:    Calculate the solid angle of an ellipse.
 *
 *      Arguments:
 *              a,b          Axis of the ellipse (x/a)^2+(y/b)^2=1
 *              x,y,z        Point at which solid angle is computed
 *      Return:
 *              Omega        solid angle
 *-----------------------------------------------------------------------*/
double SolidAngleEllipse(double a, double b, double x, double y, double z)
{
        int     Nint, maxiter, iter, j, converged;
        double  eps;
        double  omega, oldomega, theta, dtheta, weight, integrand, integral;
        double  dx, dy1, dy2, r1, r2;
        
        eps = 1e-4;            

        if (fabs(z)<eps)
        {
            if ( (x*x/a/a+y*y/b/b)<1 )
                omega = ((z>0)?(1):(-1))*(M_PI*2);
            else
                omega = 0;
            return omega;
        }
        
        omega = 0.0;
        Nint = 4;
        maxiter = 10;
        for (iter=0;iter<maxiter;iter++)
        {
            if (iter>0) oldomega = omega;
            else oldomega = 0.0;
            
            integral = 0.0;
            for (j=0;j<=Nint;j++)
            {
                if ((j==0)||(j==Nint)) weight = 0.5;
                else weight = 1.0;
                
                theta = ((j*1.0)/Nint - 0.5)*M_PI;
                dx=a*sin(theta)-x; 
                dy1=b*cos(theta)-y; dy2=-b*cos(theta)-y;
                r1=sqrt(dx*dx+dy1*dy1+z*z); r2=sqrt(dx*dx+dy2*dy2+z*z); 
                integrand = cos(theta)/(dx*dx+z*z)*( dy1/r1 - dy2/r2 );
                integral += integrand * weight;
            }
            dtheta = M_PI/Nint;
            integral *= dtheta;
            omega = a*z*integral;
            
            Nint = Nint*4;
    
            /* test convergence */
            converged = 0;
            if (iter>0)
              if (fabs(oldomega-omega) < fabs(eps*omega))
                  converged = 1;
            if (converged)
                break;
    
            if (iter==maxiter)
                fprintf(stderr,"error: solid_angle_ellipse: maxiter exceeded\n");
        }/* end of for(iter) */

        return omega;
}




