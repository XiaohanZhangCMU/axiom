/*
 * Defines python interfaces to C++ classes and methods
 */

#include <cstdio>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <pybind11/numpy.h>

#include "solver.h"

namespace bp = pybind11;
namespace axiom {

PYBIND11_MODULE(lorenz, m) {
    // optional module doc string
    m.doc() = "Lorenz solver class exposure";

    /* Expose the class Animal */
    bp::class_<Fem>(m,"Solver")
        .def(bp::init ())
        .def("fem1d_heat_steady",  &Fem::fem1d_heat_steady)
        .def("set_mesh",  &Fem::r8vec_even_new)
        .def("getU", [](Fem &m) {
            double *foo = reinterpret_cast<double*> (m.u);
	        bp::capsule free_when_done(foo, [](void *getU) {
	            double *foo = reinterpret_cast<double *>(getU);
	            std::cerr << "u[0] = " << foo[0] << "\nfreeing memory @ " << getU << std::endl;
	            delete[] foo;
	          });

	        return bp::array_t<double>(
                  {m.n, 1}, // shape
                  {1*8, 8}, // C-style contiguous strides for double
                  foo, // the data pointer
                  free_when_done); // numpy array references this parent
	        })
        .def("getX", [](Fem &m) {
            double *foo = reinterpret_cast<double*> (m.x);
	        bp::capsule free_when_done(foo, [](void *getX) {
	            double *foo = reinterpret_cast<double *>(getX);
	            std::cerr << "x[0] = " << foo[0] << "\nfreeing memory @ " << getX << std::endl;
	            delete[] foo;
	          });

	        return bp::array_t<double>(
                  {m.n, 1}, // shape
                  {1*8, 8}, // C-style contiguous strides for double
                  foo, // the data pointer
                  free_when_done); // numpy array references this parent
	        })
    ;

} /* PYBIND11_MODULE */
} /* namespace axiom */

