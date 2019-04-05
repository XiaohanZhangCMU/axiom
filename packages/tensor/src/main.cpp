#include <cstdio>
#include <iostream>
#include <pybind11/pybind11.h>

#include "tensor.hpp"
#include "axiom.hpp"

//namespace axiom { 
using namespace axiom;

int main(int argc, char *argv[]) {

    Tensor<float> t;
    //std::<<cout<<t.shape()<<std::endl;
    std::cout<<t.count()<<std::endl;
    return 0;
}
//}
