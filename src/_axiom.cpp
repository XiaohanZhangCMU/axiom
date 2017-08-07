#include <boost/python.hpp>
#include <boost/utility.hpp>

/* To prevent from deprecation warnings 
 * Need it before include arrayobject.h  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#include <cstdio>
#include <iostream>

#include "hello_zoo.hpp"
#include "hello_cuda.hpp"
#include "hello_numpy.hpp"
//#include "axiom.hpp"

/*
 * The macro Boost.Python signifies Python extension module.
 */

#define BP_REGISTER_SHARED_PTR_TO_PYTHON(PTR) do { \
  const boost::python::type_info info = \
    boost::python::type_id<shared_ptr<PTR > >(); \
  const boost::python::converter::registration* reg = \
    boost::python::converter::registry::query(info); \
  if (reg == NULL) { \
    bp::register_ptr_to_python<shared_ptr<PTR > >(); \
  } else if ((*reg).m_to_python == NULL) { \
    bp::register_ptr_to_python<shared_ptr<PTR > >(); \
  } \
} while (0)

namespace bp = boost::python;
namespace axiom {

// For Python, for now, we'll just always use float as the type.
typedef float Dtype;
const int NPY_DTYPE = NPY_FLOAT32;

// Selecting mode.
// void set_mode_cpu() { Axiom::set_mode(Axiom::CPU);  std::cerr<<"set mode to CPU"<<std::endl;}
// void set_mode_gpu() { Axiom::set_mode(Axiom::GPU);  std::cerr<<"set mode to GPU"<<std::endl;}

/* import_array() return NULL in Py3 and nothing in Py2
 * Need this wrapper for BOOST_PYTHON_MODULE if use Py3 */
void* init_numpy()
{ 
    import_array();
    return NULL;
}

void reference_contiguous_array(PyObject* in, 
                                PyArrayObject * &in_con, 
                                double* &ptr, int &count)
{
    in_con = PyArray_GETCONTIGUOUS((PyArrayObject*) in);
    ptr = (double*) PyArray_DATA(in_con);
    int num_dim = PyArray_NDIM(in_con);
    npy_intp* pdim = PyArray_DIMS(in_con);
    count = 1;
    for (int i = 0; i<num_dim; i++) count *= pdim[i];
}

void dereference(PyObject* o)
{
    Py_DECREF(o);
}

PyObject* entry_square_matrix(PyObject* input_matrix)
{
    double* ptr = NULL;
    int count = 0;
    PyArrayObject* input_contiguous_array = NULL;
    reference_contiguous_array(input_matrix, input_contiguous_array, ptr, count);
    /* create output array */
    npy_intp dst_dim[1];
    dst_dim[0] = count;
    PyObject* out_matrix = PyArray_SimpleNew(1, dst_dim, NPY_FLOAT64);
    double* ptr_out = NULL;
    PyArrayObject* output_contiguous_array = NULL;
    reference_contiguous_array(out_matrix, output_contiguous_array, ptr_out, count);
    for (int i = 0; i< count; i++)
        ptr_out[i] = ptr[i] * ptr[i];
    dereference((PyObject*)input_contiguous_array);
    dereference((PyObject*)output_contiguous_array);
    return out_matrix;
}

BOOST_PYTHON_MODULE(axiom) {

/* Expose the class Animal */
bp::class_<Animal>("Animal", bp::init<std::string const & > ())
    .def("get_address",  &Animal::get_address)
    .add_property("name", &Animal::get_name, &Animal::set_name);

/* Expose the class CudaAnimal */
bp::class_<CudaAnimal>("CudaAnimal", bp::init<std::string const & > ())
    .def("test_saxpy", &CudaAnimal::test_saxpy);

/* Expose the class NumpyAnimal */
init_numpy();
bp::class_<NumpyAnimal>("NumpyAnimal");
bp::def("square_matrix", &entry_square_matrix);

/* Expose the class Axiom */
//bp::def("set_mode_cpu", &set_mode_cpu);
//bp::def("set_mode_gpu", &set_mode_gpu);
//bp::def("set_device", &Axiom::SetDevice);
//bp::def("set_multiprocess", &Axiom::set_multiprocess);
//
//bp::class_<Tensor<Dtype>, shared_ptr<Tensor<Dtype> >, boost::noncopyable>(
//  "Tensor", bp::no_init)
//  .add_property("shape",
//      bp::make_function(
//          static_cast<const vector<int>& (Tensor<Dtype>::*)() const>(
//              &Tensor<Dtype>::shape),
//          bp::return_value_policy<bp::copy_const_reference>()))
//  .add_property("num",      &Tensor<Dtype>::num)
//  .add_property("channels", &Tensor<Dtype>::channels)
//  .add_property("height",   &Tensor<Dtype>::height)
//  .add_property("width",    &Tensor<Dtype>::width)
//  .add_property("count",    static_cast<int (Tensor<Dtype>::*)() const>(
//      &Tensor<Dtype>::count))
//  .def("reshape",           bp::raw_function(&Tensor_Reshape))
//  .add_property("data",     bp::make_function(&Tensor<Dtype>::mutable_cpu_data,
//        NdarrayCallPolicies()));
// BP_REGISTER_SHARED_PTR_TO_PYTHON(Tensor<Dtype>);
//
} /* boost_python_module */

} /* namespace axiom */

