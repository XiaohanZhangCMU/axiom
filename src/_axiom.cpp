#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <cstdio>
#include <string>
#include <vector>

#include "hello_zoo.hpp"
#include "hello_cuda.hpp"

/*
 * The macro Boost.Python signifies Python extension module.
 */

namespace bp = boost::python;
namespace _axiom {

BOOST_PYTHON_MODULE(axiom) {

// Expose the class Animal.
bp::class_<Animal>("Animal", bp::init<std::string const & > ())
	.def("get_address",  &Animal::get_address)
	.add_property("name", &Animal::get_name, &Animal::set_name);

// Expose the class CudaAnimal.
bp::class_<CudaAnimal>("CudaAnimal", bp::init<std::string const & > ())
	.def("test_saxpy",   &CudaAnimal::test_saxpy);
}

}

