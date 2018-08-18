#ifndef __TENSOR_HPP__
#define __TENSOR_HPP__

#include <algorithm>
#include <string>
#include <vector>
#include <climits>

#include "axiom.hpp"
#include "syncedmem.hpp"
#include "linalgops.hpp"

const unsigned int kMaxTensorRanks = 32; /* max tensor rank */

namespace axiom {
using namespace std;

template <typename Dtype>
void Copy(const int N, const Dtype* X, Dtype* Y) {
  if (X != Y) {
    if (Axiom::mode() == Axiom::GPU) {
#ifndef CPU_ONLY
      CUDA_CHECK(cudaMemcpy(Y, X, sizeof(Dtype) * N, cudaMemcpyDefault));
#else
      NO_GPU;
#endif
    } else {
      memcpy(Y, X, sizeof(Dtype) * N); 
    }
  }
}

/* Tensor wraps around SyncedMemory and serve as basic computation unit */
template <typename Dtype>
class Tensor {
 public:
  Tensor()
       : data_(), count_(0), capacity_(0) {}
  explicit Tensor(const vector<unsigned int>& shape);
  void Reshape(const vector<unsigned int>& shape);
  void ReshapeLike(const Tensor& other);
  inline string shape_string() const {
    ostringstream stream;
    for (unsigned int i = 0; i < shape_.size(); ++i)  stream << shape_[i] << " ";
    stream << "(" << count_ << ")";
    return stream.str();
  }

  inline const vector<unsigned int>& shape() const { return shape_; }
  inline unsigned int num_axes() const { return shape_.size(); }
  inline unsigned int count() const { return count_; }

  /* Copy from a source tensor */
  void copy(const Tensor<Dtype>& source, unsigned int reshape = 0);
  /* Assume reshaped. This is a dangerous way of assigning value to tensor */
  void reinit(const Dtype* src, unsigned int len);
  
  inline const shared_ptr<SyncedMemory>& data() const {
    CHECK(data_);
    return data_;
  }

  const Dtype* cpu_data() const;
  void set_cpu_data(Dtype* data);
  const unsigned int* gpu_shape() const;
  const Dtype* gpu_data() const;
  void set_gpu_data(Dtype* data);
  Dtype* mutable_cpu_data();
  Dtype* mutable_gpu_data();
  void Update();

  /* L1 norm */
  Dtype L1() const;
  /* L2 norm */
  Dtype L2() const;
  void ShareData(const Tensor& other);
  /* return a rank-1 tensor */
  Tensor<Dtype> operator[] (const unsigned int index) const;
  Tensor<Dtype>& operator= (const Tensor<Dtype>& otensor);
  
 protected:
  shared_ptr<SyncedMemory> data_;
  shared_ptr<SyncedMemory> shape_data_;
  vector<unsigned int> shape_;
  unsigned int count_;
  unsigned int capacity_;
//  DISABLE_COPY_AND_ASSIGN(Tensor);
  
}; /* class Tensor */
}  /* namespace axiom */
#endif /* __TENSOR_HPP__ */
