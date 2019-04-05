/*
 * Defines python interfaces to C++ classes and methods
 */

//#include <boost/python.hpp>
//#include <boost/utility.hpp>
//#include <boost/python/numpy.hpp>
#include <cstdio>
#include <iostream>
#include <memory>
#include <pybind11/pybind11.h>

#include <pybind11/operators.h>

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include "tensor.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<float>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace bp = pybind11;
namespace axiom {
typedef double Dtype;

PYBIND11_MODULE(tensor, m) {
    bp::bind_vector<std::vector<int>>(m, "VectorInt", bp::buffer_protocol());
    bp::bind_vector<std::vector<float>>(m, "VectorFloat", bp::buffer_protocol());
    bp::bind_vector<std::vector<double>>(m, "VectorDouble", bp::buffer_protocol());

    // optional module doc string
    m.doc() = "Tensor class exposure";

    /* Expose the class Animal */
    bp::class_<Tensor<Dtype> >(m,"Tensor", bp::buffer_protocol())
        .def(bp::init<> ())
        .def("Reshape", &Tensor<Dtype>::Reshape)
        .def("ReshapeLike", &Tensor<Dtype>::ReshapeLike)
        .def("shape_string", &Tensor<Dtype>::shape_string)
        .def("shape", &Tensor<Dtype>::shape)
        .def("num_axes", &Tensor<Dtype>::num_axes)
        .def("count", &Tensor<Dtype>::count)
        .def("copy",  &Tensor<Dtype>::copy)
        .def("reinit",  &Tensor<Dtype>::reinit)
        .def("data", &Tensor<Dtype>::data )
        .def("cpu_data", &Tensor<Dtype>::cpu_data)
        .def("set_cpu_data", &Tensor<Dtype>::set_cpu_data)
        .def("gpu_shape", &Tensor<Dtype>::gpu_shape)
        .def("gpu_data", &Tensor<Dtype>::gpu_data)
        .def("set_gpu_data", &Tensor<Dtype>::set_gpu_data)
        .def("L1", &Tensor<Dtype>::L1 )
        .def("L2", &Tensor<Dtype>::L2 )
        //.def("get",&Tensor<Dtype>::get)
//        .def(bp::self = bp::self)
//        .def(bp::self [] unsigned int() )

        ;

} /* PYBIND11_MODULE */
} /* namespace axiom */

