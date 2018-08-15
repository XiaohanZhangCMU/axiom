#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/python/numpy.hpp>

#include <cstdio>
#include <assert.h>
#include <iostream>

#include "hello_zoo.hpp"
#include "hello_cuda.hpp"
#include "hello_numpy.hpp"
#include "axiom.hpp"
#include "tensor.hpp"
#include "linalgops.hpp" 
#include "relax_zxcgr.hpp"

using namespace std;
using namespace axiom;

int ncalls;
double acc=1e-10;
#define Sqr(x) ((x)*(x))
void values(int n, double *x, double *f, double *g)
{
    ncalls++;
    assert(n==2);
   *f=cos(x[0]-0.5)*cos(Sqr(x[1]+1))+Sqr(x[0]-0.5)/100;
    g[0]=-sin(x[0]-0.5)*cos(Sqr(x[1]+1))+2*(x[0]-0.5)/100;
    g[1]=-cos(x[0]-0.5)*sin(Sqr(x[1]+1))*2*(x[1]+1);
//    *f=Sqr(x[0]-x[1]+0.5)+20*Sqr(x[1]+x[0]-2)+1.0;
//    g[0]=2*(x[0]-x[1]+0.5)+40*(x[1]+x[0]-2);
//    g[1]=-2*(x[0]-x[1]+0.5)+40*(x[1]+x[0]-2);
    printf("Values(%f,%f)=%f(%f,%f)\n",x[0],x[1],*f,
           g[0]/sqrt(acc),g[1]/sqrt(acc));
}

void fvalue(int *n, double *x, double *f, double *g)
{
    values(*n,x,f,g);
}

double X[2],G[2],F,W[100];

int main( )
{
#if 1
   cout<<"Hello World!"<<std::endl;
   Animal animal("I am here");
   cout<<"animal name = "<<animal.get_name()<<std::endl;
   Tensor<double> t;
   std::vector<unsigned int> shape(2);
   shape[0] = 2; shape[1] = 2;
   t.Reshape(shape);
   CudaAnimal ca("cuda animal");
   //ca.test_tensor_saxpy();
   ca.test_saxpy();
#endif


#if 1
    double x0,x1;
    int n, maxfn;
    double dfpred;

    srand(time(NULL));
    x0=rand()/100000000.0;
    x1=rand()/100000000.0;
    n=2;
    dfpred=0.02;
    maxfn=1000;

    X[0]=x0;
    X[1]=x1;
    ncalls=0;
    CGRelax(values, n, acc, maxfn, dfpred,X,G,&F, NULL);
    printf("Output:%f(%f,%f)[%f,%f] <%d>\n",
           F,X[0],X[1],G[0]/sqrt(acc),G[1]/sqrt(acc),ncalls);

    printf("\n\n\n\n\n\n");

    X[0]=x0;
    X[1]=x1;
    ncalls=0;
    CGRelax_GPU(values, n, acc, maxfn, dfpred,X,G,&F, NULL);
    printf("Output:%f(%f,%f)[%f,%f] <%d>\n",
           F,X[0],X[1],G[0]/sqrt(acc),G[1]/sqrt(acc),ncalls);

    return 0;

#endif
}
