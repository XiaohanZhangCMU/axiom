//#include <boost/thread.hpp>
#include <cmath>
#include <cstdio>
#include <ctime>
#include "axiom.hpp"

namespace axiom {

/* Make sure each thread can have different values */

//static boost::thread_specific_ptr<Axiom> thread_instance_;

Axiom& Axiom::Get() {
 // if (!thread_instance_.get()) { thread_instance_.reset(new Axiom()); }
  //return *(thread_instance_.get());
  static Axiom axiom; 
  return axiom;
}

#ifdef CPU_ONLY /* CPU-only Axiom */

Axiom::Axiom()
    : mode_(Axiom::CPU), multiprocess_(false) { }
Axiom::~Axiom() { }
void Axiom::SetDevice(const int device_id) { NO_GPU; }
void Axiom::DeviceQuery() { NO_GPU; }
bool Axiom::CheckDevice(const int device_id) { NO_GPU; return false; }
int Axiom::FindDevice(const int start_id) { NO_GPU; return -1; }

#else /* Normal GPU + CPU Axiom */

Axiom::Axiom()
    : cublas_handle_(NULL),mode_(Axiom::CPU), multiprocess_(false) {
  /* Create a cublas handler report an error if failed but
   * program still runs as one might just wanna run CPU code */
  if (cublasCreate(&cublas_handle_) != CUBLAS_STATUS_SUCCESS)
    std::cerr<<"Axiom finds no cublas_handle. cublas is not available."<<std::endl;
}

Axiom::~Axiom() {
  if (cublas_handle_) CUBLAS_CHECK(cublasDestroy(cublas_handle_));
}

void Axiom::SetDevice(const int device_id) {
  int current_device;
  CUDA_CHECK(cudaGetDevice(&current_device));
  if (current_device == device_id) {
    return;
  }
  // The call to cudaSetDevice must come before any calls to Get, which
  // may perform initialization using the GPU.
  CUDA_CHECK(cudaSetDevice(device_id));
  if (Get().cublas_handle_) CUBLAS_CHECK(cublasDestroy(Get().cublas_handle_));
  CUBLAS_CHECK(cublasCreate(&Get().cublas_handle_));
}

void Axiom::DeviceQuery() {
  cudaDeviceProp prop;
  int device;
  if (cudaSuccess != cudaGetDevice(&device)) {
    printf("No cuda device present.\n");
    return;
  }
  CUDA_CHECK(cudaGetDeviceProperties(&prop, device));
  std::cout << "Device id:                     " << device
    <<std::endl;
  std::cout << "Major revision number:         " << prop.major
    <<std::endl;
  std::cout << "Minor revision number:         " << prop.minor
    <<std::endl;
  std::cout << "Name:                          " << prop.name
    <<std::endl;
  std::cout << "Total global memory:           " << prop.totalGlobalMem
    <<std::endl;
  std::cout << "Total shared memory per block: " << prop.sharedMemPerBlock
    <<std::endl;
  std::cout << "Total registers per block:     " << prop.regsPerBlock
    <<std::endl;
  std::cout << "Warp size:                     " << prop.warpSize
    <<std::endl;
  std::cout << "Maximum memory pitch:          " << prop.memPitch
    <<std::endl;
  std::cout << "Maximum threads per block:     " << prop.maxThreadsPerBlock
    <<std::endl;
  std::cout << "Maximum dimension of block:    "
      << prop.maxThreadsDim[0] << ", " << prop.maxThreadsDim[1] << ", "
      << prop.maxThreadsDim[2]<<std::endl;
  std::cout << "Maximum dimension of grid:     "
      << prop.maxGridSize[0] << ", " << prop.maxGridSize[1] << ", "
      << prop.maxGridSize[2]<<std::endl;
  std::cout << "Clock rate:                    " << prop.clockRate
    <<std::endl;
  std::cout << "Total constant memory:         " << prop.totalConstMem
    <<std::endl;
  std::cout << "Texture alignment:             " << prop.textureAlignment
    <<std::endl;
  std::cout << "Concurrent copy and execution: "
      << (prop.deviceOverlap ? "Yes" : "No")<<std::endl;
  std::cout << "Number of multiprocessors:     " << prop.multiProcessorCount
    <<std::endl;
  std::cout << "Kernel execution timeout:      "
      << (prop.kernelExecTimeoutEnabled ? "Yes" : "No")<<std::endl;
  return;
}

bool Axiom::CheckDevice(const int device_id) {
  // This function checks the availability of GPU #device_id.
  // It attempts to create a context on the device by calling cudaFree(0).
  // cudaSetDevice() alone is not sufficient to check the availability.
  // It lazily records device_id, however, does not initialize a
  // context. So it does not know if the host thread has the permission to use
  // the device or not.
  //
  // In a shared environment where the devices are set to EXCLUSIVE_PROCESS
  // or EXCLUSIVE_THREAD mode, cudaSetDevice() returns cudaSuccess
  // even if the device is exclusively occupied by another process or thread.
  // Cuda operations that initialize the context are needed to check
  // the permission. cudaFree(0) is one of those with no side effect,
  // except the context initialization.
  bool r = ((cudaSuccess == cudaSetDevice(device_id)) &&
            (cudaSuccess == cudaFree(0)));
  // reset any error that may have occurred.
  cudaGetLastError();
  return r;
}

int Axiom::FindDevice(const int start_id) {
  // This function finds the first available device by checking devices with
  // ordinal from start_id to the highest available value. In the
  // EXCLUSIVE_PROCESS or EXCLUSIVE_THREAD mode, if it succeeds, it also
  // claims the device due to the initialization of the context.
  int count = 0;
  CUDA_CHECK(cudaGetDeviceCount(&count));
  for (int i = start_id; i < count; i++) {
    if (CheckDevice(i)) return i;
  }
  return -1;
}
#endif /* CPU + GPU */
} /* namespace axiom */
