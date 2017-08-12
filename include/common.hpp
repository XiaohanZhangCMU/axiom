#ifndef __COMMON_H__
#define __COMMON_H__

#include <boost/shared_ptr.hpp> 
#include <boost/thread.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

//#define CPU_ONLY

#define LOG(...) std::cout , __VA_ARGS__ , std::endl

#define CHECK_EQ(x, y) \
  do { \
    if(x != y){ \
     printf("CHECK_EQ failure %s:%d\n",__FILE__,__LINE__); \
     exit(1);} \
  } while (0)

#define CHECK_LE(x, y) \
  do { \
    if(x > y){ \
     printf("CHECK_LE failure %s:%d\n",__FILE__,__LINE__); \
     exit(1);} \
  } while (0)
#define CHECK_LT(x, y) \
  do { \
    if(x >= y){ \
     printf("CHECK_LE failure %s:%d\n",__FILE__,__LINE__); \
     exit(1);} \
  } while (0)
#define CHECK_GE(x, y) \
  do { \
    if(x < y){ \
     printf("CHECK_GE failure %s:%d\n",__FILE__,__LINE__); \
     exit(1);} \
  } while (0)
#define CHECK_GT(x, y) \
  do { \
    if(x <= y){ \
     printf("CHECK_GE failure %s:%d\n",__FILE__,__LINE__); \
     exit(1);} \
  } while (0)

#define CHECK(x) \
  do { \
    if( !(x) ){ \
     printf("failure %s:%d\n",__FILE__,__LINE__); \
     exit(1);} \
  } while (0)

/* Disable the copy and assignment operator for a class */
#define DISABLE_COPY_AND_ASSIGN(classname) \
private:\
  classname(const classname&);\
  classname& operator=(const classname&)

/* Instantiate a class with float and double specifications */
#define INSTANTIATE_CLASS(classname) \
  char gInstantiationGuard##classname; \
  template class classname<float>; \
  template class classname<double>

#define NOT_IMPLEMENTED printf("NotImplemented. %s:%d\n", __FILE__,__LINE__);

#ifdef CPU_ONLY  /* CPU-only Axiom */

#define NO_GPU printf("Cannot use GPU in CPU-only Axiom: check mode");

#else  /* Normal GPU + CPU Axiom */

#include <cublas_v2.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <driver_types.h>  

/* CUDA: various checks for different function calls */
#define CUDA_CHECK(condition) \
  /* Code block avoids redefinition of cudaError_t error */ \
  do { \
    cudaError_t error = condition; \
    if(error != cudaSuccess) \
     printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(error)); \
  } while (0)

//     printf("cublas failure %s:%d: '%s'\n",__FILE__,__LINE__,axiom::CublasGetErrorString(status)); 

#define CUBLAS_CHECK(condition) \
  do { \
    cublasStatus_t status = condition; \
    if(status != CUBLAS_STATUS_SUCCESS) \
     printf("cublas failure %s:%d: '%s'\n",__FILE__,__LINE__,"CUBLAS_STATUS_SUCCESS"); \
  } while (0)

/* CUDA: grid stride looping */
#define CUDA_KERNEL_LOOP(i, n) \
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; \
       i < (n); \
       i += blockDim.x * gridDim.x)

/* CUDA: check for error after kernel execution and exit loudly if there is one */
#define CUDA_POST_KERNEL_CHECK CUDA_CHECK(cudaPeekAtLastError())

namespace axiom {

/* CUDA: library error reporting */
//const char* CublasGetErrorString(cublasStatus_t error) {
//    switch (error) {
//    case CUBLAS_STATUS_SUCCESS:
//        return "CUBLAS_STATUS_SUCCESS";
//    case CUBLAS_STATUS_NOT_INITIALIZED:
//        return "CUBLAS_STATUS_NOT_INITIALIZED";
//    case CUBLAS_STATUS_ALLOC_FAILED:
//        return "CUBLAS_STATUS_ALLOC_FAILED";
//    case CUBLAS_STATUS_INVALID_VALUE:
//        return "CUBLAS_STATUS_INVALID_VALUE";
//    case CUBLAS_STATUS_ARCH_MISMATCH:
//        return "CUBLAS_STATUS_ARCH_MISMATCH";
//    case CUBLAS_STATUS_MAPPING_ERROR:
//        return "CUBLAS_STATUS_MAPPING_ERROR";
//    case CUBLAS_STATUS_EXECUTION_FAILED:
//        return "CUBLAS_STATUS_EXECUTION_FAILED";
//    case CUBLAS_STATUS_INTERNAL_ERROR:
//        return "CUBLAS_STATUS_INTERNAL_ERROR";
//#if CUDA_VERSION >= 6000
//    case CUBLAS_STATUS_NOT_SUPPORTED:
//        return "CUBLAS_STATUS_NOT_SUPPORTED";
//#endif
//#if CUDA_VERSION >= 6050
//    case CUBLAS_STATUS_LICENSE_ERROR:
//        return "CUBLAS_STATUS_LICENSE_ERROR";
//#endif
//    }
//    return "Unknown cublas status";
//}

/* CUDA: use 512 threads per block */
const int AXIOM_CUDA_NUM_THREADS = 512; 
/* CUDA: number of blocks for threads */
inline int AXIOM_GET_BLOCKS(const int N) {
  return (N + AXIOM_CUDA_NUM_THREADS - 1) / AXIOM_CUDA_NUM_THREADS;
}
} /* namespace axiom */ 
#endif /* GPU+CPU */
#endif /* __COMMON_H__ */
