#ifndef __HELLO_ZOO_HPP__
#define __HELLO_ZOO_HPP__

#include <cstdint>
#include <string>
#include <vector>

namespace axiom {

class Animal {

public:	
    /* Animal constructors */	
    Animal(std::string const & in_name) : m_name(in_name) {}	
    Animal(Animal const & in_other) : m_name(in_other.m_name) {}
    Animal & operator = (Animal const & in_other) {
    	this->m_name = in_other.m_name;
    	return *this;
    }

    /* C++ function to expose to python */
    uintptr_t get_address() const; 
    /* Getter of the name property */
    std::string get_name() const;
    /* Setter of the name property */
    void set_name(std::string const & in_name);
private:	
    std::string m_name; /* some property */

}; /* class Animal */ 
} /* namespace axiom */
#endif /* __HELLO_ZOO_HPP__ */
