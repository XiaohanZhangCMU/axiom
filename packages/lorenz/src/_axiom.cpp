/*
 * Defines python interfaces to C++ classes and methods
 */

#include <cstdio>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <pybind11/numpy.h>

#include "sim.h"

namespace bp = pybind11;
namespace axiom {

PYBIND11_MODULE(lorenz, m) {
    // optional module doc string
    m.doc() = "Sim class exposure";

    /* Expose the class Animal */
    bp::class_<Sim>(m,"Sim")
        .def(bp::init ())
    //    .def("run", &Sim::run)
    ;
    m.def("fillVBO", fillVBO);

} /* PYBIND11_MODULE */
} /* namespace axiom */

