#ifndef __FEM1D_HPP__
#define __FEM1D_HPP__

#include <cstdint>
#include <string>
#include <vector>

namespace axiom {

class Fem {

public:	
    /* constructors */	
    Fem() {}	

    /* C++ function to expose to python */
    void fem1d_heat_steady ( int n, double a, double b, double ua, double ub);

    void r8vec_even_new ( int n, double alo, double ahi );

    int n;
    double* u;
    double* x; //mesh

private:
    /* Private functions */
    int i4_max ( int i1, int i2 );
    int i4_min ( int i1, int i2 );
    int *i4vec_zero_new ( int n );
    double k ( double x );
    double f ( double x );
    double r8_abs ( double x );
    void r8mat_print ( int m, int n, double a[], std::string title );
    void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
      int jhi, std::string title );
    double *r8mat_solve2 ( int n, double a[], double b[], int *ierror );
    double *r8mat_zero_new ( int m, int n );
    double *r8vec_zero_new ( int n );
    void timestamp ( );


}; /* class Fem */ 
} /* namespace axiom */
#endif /* __FEM1D_HPP__ */
