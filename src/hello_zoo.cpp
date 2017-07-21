#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <cstdint>
#include <string>
#include <vector>

/*
 * C++ function to expose to python
 */
const std::string hello() {
	return std::string("hello, zoo");
}

/*
 * C++ class represents animals in zoo
 */
class Animal {

public:
	Animal(std::string const & in_name) : m_name(in_name) {}	

	Animal(Animal const & in_other) : m_name(in_other.m_name) {}

	Animal & operator = (Animal const & in_other) {
		this->m_name = in_other.m_name;
		return *this;
	}

	// Utility method to get address of the instance
	uintptr_t get_address() const {
		return reinterpret_cast<uintptr_t>(this);
	}

	// Getter of the name property
	std::string get_name() const {
		return this->m_name;
	}

	// Setter of the name property
	void set_name(std::string const & in_name) {
		this->m_name = in_name;
	}
private:	
	//The property
	std::string m_name;
};

/*
 * The macro Boost.Python signifies Python extension module.
 */
BOOST_PYTHON_MODULE(axiom) {
	using namespace boost::python;

	// Expose the function hello().
	def("hello", hello);

	// Expose the class Animal.
	class_<Animal>("Animal", 
			init<std::string const & > ())
		.def("get_address", &Animal::get_address)
		.add_property("name", &Animal::get_name, &Animal::set_name);
}
