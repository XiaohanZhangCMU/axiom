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
        .def("run_fine_ODE", &Sim::run_fine_ODE)
        .def("trajectory", [](Sim &m) {
            double *foo = reinterpret_cast<double*> (m.trajectory);
	        bp::capsule free_when_done(foo, [](void *trajectory) {
	            //double *foo = reinterpret_cast<double *>(trajectory);
	            //std::cerr << "SR [0] = " << foo[0] << "\nfreeing memory @ " << trajectory << std::endl;
	            //delete[] foo;
	          });

	        return bp::array_t<double>(
                  {m.points, 3}, // shape
                  {3*8, 8}, // C-style contiguous strides for double
                  foo, // the data pointer
                  free_when_done); // numpy array references this parent
	        })
        .def_readwrite("points",&Sim::points)
        .def_readwrite("t",&Sim::t)
        .def_readwrite("dt",&Sim::dt)
    ;

} /* PYBIND11_MODULE */
} /* namespace axiom */

