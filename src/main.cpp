#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/python/numpy.hpp>

#include <cstdio>
#include <iostream>

#include "hello_zoo.hpp"
#include "hello_cuda.hpp"
#include "hello_numpy.hpp"
#include "axiom.hpp"
#include "tensor.hpp"
#include "linalgops.hpp" 

using namespace std;
using namespace axiom;

int main( )
{
   cout<<"Hello World!"<<std::endl;
   Animal animal("I am here");
   cout<<"animal name = "<<animal.get_name()<<std::endl;
   Tensor<double> t;
   std::vector<unsigned int> shape(2);
   shape[0] = 2; shape[1] = 2;
   t.Reshape(shape);
}
