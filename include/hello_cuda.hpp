#ifndef __HELLO_CUDA_HPP__
#define __HELLO_CUDA_HPP__

#include <stdio.h>
#include <cstdint>
#include <string>
#include <vector>

namespace _axiom {

/*
 * Adapted from https://gist.github.com/kazimuth/f9952810d8117f782be71279c0f16f6e
 */

class CudaAnimal {

public:	
	/*
	 * Animal constructors
	 */	
	CudaAnimal(std::string const & in_name) : m_name(in_name) {}	
	CudaAnimal(CudaAnimal const & in_other) : m_name(in_other.m_name) {}
	CudaAnimal & operator = (CudaAnimal const & in_other) {
		this->m_name = in_other.m_name;
		return *this;
	}

	int test_saxpy(void);

private:	
	//The property
	std::string m_name;
}; 

}
#endif
