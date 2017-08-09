#ifndef __SYNCEDMEM_HPP__
#define __SYNCEDMEM_HPP__

#include <cstdlib>
#include "common.hpp"
#include "axiom.hpp"

namespace axiom {

inline void Malloc(void** ptr, size_t size, bool* use_cuda) {
#ifndef CPU_ONLY
  if (Axiom::mode() == Axiom::GPU) {
    CUDA_CHECK(cudaMallocHost(ptr, size));
    *use_cuda = true;
    return;
  }
#endif
  *ptr = malloc(size);
  *use_cuda = false;
  CHECK(ptr);
}

inline void Free(void* ptr, bool use_cuda) {
#ifndef CPU_ONLY
  if (use_cuda) {
    CUDA_CHECK(cudaFreeHost(ptr));
    return;
  }
#endif
  free(ptr);
}

inline void HostMemset(const size_t N, const int alpha, void* X) {
  memset(X, alpha, N);
}
inline void DeviceMemset(const size_t N, const int alpha, void* X) {
#ifndef CPU_ONLY
  CUDA_CHECK(cudaMemset(X, alpha, N));
#else
  NO_GPU;
#endif
}

inline void DeviceMemcpy(const size_t N, const void* X, void* Y) {
#ifndef CPU_ONLY
  if (X != Y) {
    CUDA_CHECK(cudaMemcpy(Y, X, N, cudaMemcpyDefault));
  }
#else
  NO_GPU
#endif
}

/* Memory synchronization between Host (CPU) and Device (GPU) */

class SyncedMemory {
public:
  SyncedMemory();
  explicit SyncedMemory(size_t size);
  ~SyncedMemory();
  const void* cpu_data();
  void set_cpu_data(void* data);
  const void* gpu_data();
  void set_gpu_data(void* data);
  void* mutable_cpu_data();
  void* mutable_gpu_data();
  enum SyncedHead { UNINITIALIZED, HEAD_AT_CPU, HEAD_AT_GPU, SYNCED };
  SyncedHead head() { return head_; }
  size_t size() { return size_; }

#ifndef CPU_ONLY
  void async_gpu_push(const cudaStream_t& stream);
#endif

private:
  void check_device();
  void to_cpu();
  void to_gpu();
  void* cpu_ptr_;
  void* gpu_ptr_;
  size_t size_;
  SyncedHead head_;
  bool own_cpu_data_;
  bool cpu_malloc_use_cuda_;
  bool own_gpu_data_;
  int device_;
  DISABLE_COPY_AND_ASSIGN(SyncedMemory);

}; /* class SyncedMemory */
}  /* namespace axiom */
#endif /* __SYNCEDMEM_HPP__ */
