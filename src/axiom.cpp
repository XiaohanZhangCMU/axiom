#include <boost/thread.hpp>
#include <cmath>
#include <cstdio>
#include <ctime>
#include "axiom.hpp"

namespace axiom {

/* Make sure each thread can have different values */
static boost::thread_specific_ptr<Axiom> thread_instance_;

Axiom& Axiom::Get() {
  if (!thread_instance_.get()) {
    thread_instance_.reset(new Axiom());
  }
  return *(thread_instance_.get());
}


#ifdef CPU_ONLY  // CPU-only Axiom.

Axiom::Axiom()
    : mode_(Axiom::CPU), multiprocess_(false) { }

Axiom::~Axiom() { }

#else  // Normal GPU + CPU Axiom.

Axiom::Axiom()
    : cublas_handle_(NULL),
    mode_(Axiom::CPU), multiprocess_(false) {
  /* Create a cublas handler report an error if failed but
   * program still runs as one might just wanna run CPU code */
  if (cublasCreate(&cublas_handle_) != CUBLAS_STATUS_SUCCESS) {
    std::cerr<<"Cannot find Cublas handle.";
    std::cerr<<"Cublas is not available..."<<std::endl;
  }
}

Axiom::~Axiom() {
  if (cublas_handle_) CUBLAS_CHECK(cublasDestroy(cublas_handle_));
}

#endif /* CPU + GPU */
} /* namespace axiom */
