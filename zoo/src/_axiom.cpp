/*
 * Defines python interfaces to C++ classes and methods
 */

#include <cstdio>
#include <iostream>
#include <pybind11/pybind11.h>

#include "hello_zoo.hpp"

namespace bp = pybind11;
namespace axiom {

PYBIND11_MODULE(axiomlib, m) {
    // optional module doc string
    m.doc() = "Animal class exposure";

    /* Expose the class Animal */
    bp::class_<Animal>(m,"Animal")
        .def(bp::init<std::string const & > ())
        .def("get_address",  &Animal::get_address)
        .def("get_name",     &Animal::get_name)
        ;

} /* PYBIND11_MODULE */
} /* namespace axiom */

