/*
 * Defines python interfaces to C++ classes and methods
 */

#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/python/numpy.hpp>
#include <cstdio>
#include <iostream>
#include <pybind11/pybind11.h>

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include "tensor.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<int>);

namespace bp = pybind11;
namespace axiom {
typedef double Dtype;

PYBIND11_MODULE(tensor, m) {
    bp::bind_vector<std::vector<int>>(m, "VectorInt", bp::buffer_protocol());
    // optional module doc string
    m.doc() = "Tensor class exposure";

    /* Expose the class Animal */
    bp::class_<Tensor<Dtype> >(m,"Tensor")
        .def(bp::init<> ())
        //.def("shape", &Tensor<Dtype>::shape<std::vector<int>>, bp::return_value_policy::reference)
        .def("shape", &Tensor<Dtype>::shape, bp::return_value_policy::reference)
        ;

} /* PYBIND11_MODULE */
} /* namespace axiom */

