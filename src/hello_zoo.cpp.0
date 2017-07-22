#include <hello_zoo.hpp>

/*
 * C++ function to expose to python
 */

namespace _axiom {

// Utility method to get address of the instance
uintptr_t Animal::get_address() const {
	return reinterpret_cast<uintptr_t>(this);
}

// Getter of the name property
std::string Animal::get_name() const {
	return this->m_name;
}

// Setter of the name property
void Animal::set_name(std::string const & in_name) {
	this->m_name = in_name;
}

}
