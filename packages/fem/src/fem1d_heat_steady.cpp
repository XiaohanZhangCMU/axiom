# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

#include "fem1d_heat_steady.hpp"

namespace axiom {

//****************************************************************************80

void Fem::fem1d_heat_steady ( int n, double a, double b, double ua, double ub)

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_HEAT_STEADY solves the steady 1D heat equation with finite elements.
//
//  Discussion:
//
//    The program uses the finite element method, with piecewise linear basis
//    functions to solve the steady state heat equation in one dimension.
//
//    The problem is defined on the region A <= x <= B.
//
//    The following differential equation is imposed between A and B:
//
//      - d/dx k(x) du/dx = f(x)
//
//    where k(x) and f(x) are given functions.
//
//    At the boundaries, the following conditions are applied:
//
//      u(A) = UA
//      u(B) = UB
//
//    A set of N equally spaced nodes is defined on this
//    interval, with A = X(1) < X(2) < ... < X(N) = B.
//
//    At each node I, we associate a piecewise linear basis function V(I,X),
//    which is 0 at all nodes except node I.  This implies that V(I,X) is
//    everywhere 0 except that
//
//    for X(I-1) <= X <= X(I):
//
//      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) ) 
//
//    for X(I) <= X <= X(I+1):
//
//      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )
//
//    We now assume that the solution U(X) can be written as a linear
//    sum of these basis functions:
//
//      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)
//
//    where U(X) on the left is the function of X, but on the right,
//    is meant to indicate the coefficients of the basis functions.
//
//    To determine the coefficient U(J), we multiply the original
//    differential equation by the basis function V(J,X), and use
//    integration by parts, to arrive at the I-th finite element equation:
//
//        Integral K(X) * U'(X) * V'(I,X) dx = Integral F(X) * V(I,X) dx
//
//    We note that the functions U(X) and U'(X) can be replaced by
//    the finite element form involving the linear sum of basis functions,
//    but we also note that the resulting integrand will only be nonzero
//    for terms where J = I - 1, I, or I + 1.
//
//    By writing this equation for basis functions I = 2 through N - 1,
//    and using the boundary conditions, we have N linear equations
//    for the N unknown coefficients U(1) through U(N), which can
//    be easily solved.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double A, B, the left and right endpoints.
//
//    Input, double UA, UB, the prescribed value of U at A and B.
//
//    Input, double K ( double X ), evaluates k(x);
//
//    Input, double F ( double X ), evaluates f(x);
//
//    Input, double X[N], the mesh points.
//
//    Output, double FEM1D_HEAT_STEADY[N], the finite element coefficients, 
//    which are also the value of the computed solution at the mesh points.
//
{
# define QUAD_NUM 2

  int i, ierror, q;
  int quad_num = QUAD_NUM;
  double abscissa[QUAD_NUM] = {
    -0.577350269189625764509148780502,
    +0.577350269189625764509148780502 };
  double weight[QUAD_NUM] = { 1.0, 1.0 };
  double al, am, ar, wq, vlp, vm, vmp, vrp, xl, xm,xq, xr,bm, fxq, kxq; 
  double *amat, *bvec;

  this->n = n;
  amat = r8mat_zero_new ( n, n );
  bvec = r8vec_zero_new ( n );
//
//  Equation 1 is the left boundary condition, U(A) = UA;
//
  amat[0+0*n] = 1.0;
  bvec[0] = ua;

  for ( i = 1; i < n - 1; i++ ) {
//  Get the left, right and middle coordinates.
    xl = x[i-1];
    xm = x[i];
    xr = x[i+1];

//  Make temporary variables for A(I,I-1), A(I,I), A(I,I+1) and B(I).
    al = am = ar = bm = 0.0;

//  We approximate the integrals by using a weighted sum of
//  the integrand values at quadrature points.
    for ( q = 0; q < quad_num; q++ )
    {
//
//  Integrate over the LEFT interval, between XL and XM, where:
//
//  VL(X) = ( XM - X       ) / ( XM - XL )
//  VM(X) = (      X  - XL ) / ( XM - XL )
//  VR(X) = 0
//
//  VL'(X) =             - 1 / ( XM - XL )
//  VM'(X) =             + 1 / ( XM - XL ) 
//  VR'(X) = 0
//
      xq = ( ( 1.0 - abscissa[q] ) * xl 
           + ( 1.0 + abscissa[q] ) * xm ) 
           /   2.0;

      wq = weight[q] * ( xm - xl ) / 2.0;

      //vl =  ( xm - xq ) / ( xm - xl );
      vlp =      - 1.0  / ( xm - xl );

      vm =  ( xq - xl ) / ( xm - xl );
      vmp =      + 1.0  / ( xm - xl );

      //vr =  0.0;
      vrp = 0.0;

      kxq = k ( xq );
      fxq = f ( xq );

      al = al + wq * ( kxq * vlp * vmp );
      am = am + wq * ( kxq * vmp * vmp );
      ar = ar + wq * ( kxq * vrp * vmp );
      bm = bm + wq * ( fxq * vm );

//  Integrate over the RIGHT interval, between XM and XR, where:
//  VL(X) = 0
//  VM(X) = ( XR - X       ) / ( XR - XM )
//  VR(X) = (      X  - XM ) / ( XR - XM )
//
//  VL'(X) = 0
//  VM'(X) =             - 1 / ( XR - XM )
//  VR'(X) =             + 1 / ( XR - XM ) 
      xq = ( ( 1.0 - abscissa[q] ) * xm 
           + ( 1.0 + abscissa[q] ) * xr ) 
           /   2.0;

      wq = weight[q] * ( xr - xm ) / 2.0;

      //vl = 0.0;
      vlp = 0.0;

      vm = ( xr - xq ) / ( xr - xm );
      vmp =     - 1.0  / ( xr - xm );

      //vr = ( xq - xm ) / ( xr - xm );
      vrp =      1.0   / ( xr - xm );

      kxq = k ( xq );
      fxq = f ( xq );

      al = al + wq * ( kxq * vlp * vmp );
      am = am + wq * ( kxq * vmp * vmp );
      ar = ar + wq * ( kxq * vrp * vmp );
      bm = bm + wq * ( fxq * vm );
    }
    amat[i+(i-1)*n] = al;
    amat[i+ i   *n] = am;
    amat[i+(i+1)*n] = ar;

    bvec[i] = bm;
  }
//  Equation N is the right boundary condition, U(B) = UB;
  amat[n-1+(n-1)*n] = 1.0;
  bvec[n-1] = ub;

//  Solve the linear system.
  u = r8mat_solve2 ( n, amat, bvec, &ierror );

  delete [] amat;
  delete [] bvec;

//  return u;
# undef QUAD_NUM
}

void Fem::r8vec_even_new ( int n, double alo, double ahi ) {
  this->x = new double[n];
  if ( n == 1 ) x[0] = 0.5 * ( alo + ahi );
  else { 
    for ( int i = 0; i < n; i++ )
      this->x[i] = ( ( double ) ( n - i - 1 ) * alo 
             + ( double ) (     i     ) * ahi ) 
             / ( double ) ( n     - 1 );
  }
}

double Fem::k ( double x ) {
  double value = 1.0;
  return value;
}
double Fem::f ( double x ) {
  double value = x * ( x + 3.0 ) * exp ( x );
  return value;
}

double * Fem::r8mat_solve2 ( int n, double a[], double b[], int *ierror )

//****************************************************************************80
//
//  Solves A*X = B
//
//  Parameters:
//
//    Input, int N, the number of equations.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix to be inverted.
//    On output, A has been overwritten.
//
//    Input/output, double B[N].
//    On input, B is the right hand side of the system.
//    On output, B has been overwritten.
//
//    Output, int *IERROR.
//    0, no error detected.
//    1, consistent singularity.
//    2, inconsistent singularity.
//
//    Output, double R8MAT_SOLVE2[N], the solution of the linear system.
//
{
  double amax, *x;
  int i, imax, j, k, *piv;

  *ierror = 0;

  piv = i4vec_zero_new ( n );
  x = r8vec_zero_new ( n );
//
//  Process the matrix.
//
  for ( k = 1; k <= n; k++ )
  {
//
//  In column K:
//    Seek the row IMAX with the properties that:
//      IMAX has not already been used as a pivot;
//      A(IMAX,K) is larger in magnitude than any other candidate.
//
    amax = 0.0;
    imax = 0;
    for ( i = 1; i <= n; i++ )
    {
      if ( piv[i-1] == 0 )
      {
        if ( amax < r8_abs ( a[i-1+(k-1)*n] ) )
        {
          imax = i;
          amax = r8_abs ( a[i-1+(k-1)*n] );
        }
      }
    }
//
//  If you found a pivot row IMAX, then,
//    eliminate the K-th entry in all rows that have not been used for pivoting.
//
    if ( imax != 0 )
    {
      piv[imax-1] = k;
      for ( j = k+1; j <= n; j++ )
      {
        a[imax-1+(j-1)*n] = a[imax-1+(j-1)*n] / a[imax-1+(k-1)*n];
      }
      b[imax-1] = b[imax-1] / a[imax-1+(k-1)*n];
      a[imax-1+(k-1)*n] = 1.0;

      for ( i = 1; i <= n; i++ )
      {
        if ( piv[i-1] == 0 )
        {
          for ( j = k+1; j <= n; j++ )
          {
            a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - a[i-1+(k-1)*n] * a[imax-1+(j-1)*n];
          }
          b[i-1] = b[i-1] - a[i-1+(k-1)*n] * b[imax-1];
          a[i-1+(k-1)*n] = 0.0;
        }
      }
    }
  }
//
//  Now, every row with nonzero IPIV begins with a 1, and
//  all other rows are all zero.  Begin solution.
//
  for ( j = n; 1 <= j; j-- )
  {
    imax = 0;
    for ( k = 1; k <= n; k++ )
    {
      if ( piv[k-1] == j )
      {
        imax = k;
      }
    }

    if ( imax == 0 )
    {
      x[j-1] = 0.0;

      if ( b[j-1] == 0.0 )
      {
        *ierror = 1;
        cout << "\n";
        cout << "R8MAT_SOLVE2 - Warning:\n";
        cout << "  Consistent singularity, equation = " << j << "\n";
      }
      else
      {
        *ierror = 2;
        cout << "\n";
        cout << "R8MAT_SOLVE2 - Warning:\n";
        cout << "  Inconsistent singularity, equation = " << j << "\n";
      }
    }
    else
    {
      x[j-1] = b[imax-1];

      for ( i = 1; i <= n; i++ )
      {
        if ( i != imax )
        {
          b[i-1] = b[i-1] - a[i-1+(j-1)*n] * x[j-1];
        }
      }
    }
  }

  delete [] piv;

  return x;
}

//
// ********************** Helper Functions ***************************
//

int Fem::i4_max ( int i1, int i2 ) { return i2<i1?i1:i2; }
int Fem::i4_min ( int i1, int i2 ) { return i1<i2?i1:i2; }
double Fem::r8_abs ( double x ) { return abs(x); }
int *Fem::i4vec_zero_new ( int n ) {
  int *a = new int[n];
  for (int i = 0; i < n; i++ ) a[i] = 0;
  return a;
}
double *Fem::r8mat_zero_new ( int m, int n ) {
  double *a = new double[m*n];
  for (int j = 0; j < n; j++ ) for (int i = 0; i < m; i++ ) a[i+j*m] = 0.0;
  return a;
}
double *Fem::r8vec_zero_new ( int n ) {
  double *a = new double[n];
  for (int i = 0; i < n; i++ ) a[i] = 0.0;
  return a;
}

//
// ********************** Printing Utilities ***************************
//
void Fem::r8mat_print ( int m, int n, double a[], std::string title ) {
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );
}
void Fem::r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, std::string title ) {
# define INCX 5

  int i, i2hi, i2lo, j, j2hi, j2lo;
  cout << "\n"<< title << "\n";

  if ( m <= 0 || n <= 0 ) {
    cout << "\n" << "  (None)\n";
    return;
  }

//  Print the columns of the matrix, in strips of 5.
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX ) {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );
    cout << "\n";

//  For each column J in the current range, write the header.
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
      cout << setw(7) << j - 1 << "       ";
    cout << "\n"<< "  Row\n"<< "\n";

//  Determine the range of the rows in this strip.
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ ) {
//  Print out (up to) 5 entries in row I, that lie in the current strip.
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ ) 
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      cout << "\n";
    }
  }
  return;
# undef INCX
}

void Fem::timestamp ( ) { // prints current time stamp
# define TIME_SIZE 40
  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  //size_t len;
  std::time_t now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );
  std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );
  std::cout << time_buffer << "\n";
  return;
# undef TIME_SIZE
}

}
