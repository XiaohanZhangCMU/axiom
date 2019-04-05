#include "tensor.hpp"

namespace axiom {

template <typename Dtype>
void Tensor<Dtype>::Reshape(const vector<unsigned int>& shape) {
  CHECK_LE(shape.size(), kMaxTensorRanks);
  count_ = 1;
  shape_.resize(shape.size());
  if (!shape_data_ || shape_data_->size() < shape.size() * sizeof(unsigned int))
    shape_data_.reset(new SyncedMemory(shape.size() * sizeof(unsigned int)));

  unsigned int* shape_data = static_cast<unsigned int*>(shape_data_->mutable_cpu_data());
  for (unsigned int i = 0; i < shape.size(); ++i) {
    CHECK_GE(shape[i], 0);
    /* tensor size cannot exceeds INT_MAX */
    if (count_ != 0) CHECK_LE(shape[i], INT_MAX / count_) ;
    count_ *= shape[i];
    shape_[i] = shape[i];
    shape_data[i] = shape[i];
  }
  if (count_ > capacity_) {
    capacity_ = count_;
    data_.reset(new SyncedMemory(capacity_ * sizeof(Dtype)));
  }
}

template <typename Dtype>
void Tensor<Dtype>::ReshapeLike(const Tensor<Dtype>& other) {
  Reshape(other.shape());
}

template <typename Dtype>
Tensor<Dtype>::Tensor(const vector<unsigned int>& shape)
  /* capacity_ must be initialized before calling Reshape */
  : capacity_(0) {
  Reshape(shape);
}

template <typename Dtype>
const unsigned int* Tensor<Dtype>::gpu_shape() const {
  CHECK(shape_data_);
  return (const unsigned int*)shape_data_->gpu_data();
}

template <typename Dtype>
const Dtype* Tensor<Dtype>::cpu_data() const {
  CHECK(data_);
  return (const Dtype*)data_->cpu_data();
}

/* shallow copy on cpu */
template <typename Dtype>
void Tensor<Dtype>::set_cpu_data(Dtype* data) {
  CHECK(data);
  /* CPU and GPU sizes must be equal */
  size_t size = count_ * sizeof(Dtype);
  if (data_->size() != size) {
    data_.reset(new SyncedMemory(size));
  }
  data_->set_cpu_data(data);
}

/* return data on gpu */
template <typename Dtype>
const Dtype* Tensor<Dtype>::gpu_data() const {
  CHECK(data_);
  return (const Dtype*)data_->gpu_data();
}

/* shallow copy on gpu */
template <typename Dtype>
void Tensor<Dtype>::set_gpu_data(Dtype* data) {
  CHECK(data);
  /* CPU and GPU sizes must be equal */
  size_t size = count_ * sizeof(Dtype);
  if (data_->size() != size) {
    data_.reset(new SyncedMemory(size));
  }
  data_->set_gpu_data(data);
}

template <typename Dtype>
Dtype* Tensor<Dtype>::mutable_cpu_data() {
  CHECK(data_);
  return static_cast<Dtype*>(data_->mutable_cpu_data());
}

template <typename Dtype>
Dtype* Tensor<Dtype>::mutable_gpu_data() {
  CHECK(data_);
  return static_cast<Dtype*>(data_->mutable_gpu_data());
}

template <typename Dtype>
void Tensor<Dtype>::ShareData(const Tensor& other) {
  CHECK_EQ(count_, other.count());
  data_ = other.data();
}

template <> unsigned int Tensor<unsigned int>::L1() const {
  NOT_IMPLEMENTED;
  return 0;
}

template <> int Tensor<int>::L1() const {
  NOT_IMPLEMENTED;
  return 0;
}

template <typename Dtype>
Dtype Tensor<Dtype>::L1() const {
  if (!data_) { return 0; }
  switch (data_->head()) {
  case SyncedMemory::HEAD_AT_CPU:
    return cpu_asum(count_, cpu_data());
  case SyncedMemory::HEAD_AT_GPU:
  case SyncedMemory::SYNCED:
#ifndef CPU_ONLY
  {
    Dtype asum;
    gpu_asum(count_, gpu_data(), &asum);
    return asum;
  }
#else
    NO_GPU;
#endif
  case SyncedMemory::UNINITIALIZED:
    return 0;
  default:
    std::cerr << "Unknown SyncedMemory head state: " << data_->head()<<std::endl;
  }
  return 0;
}

template <> unsigned int Tensor<unsigned int>::L2() const {
  NOT_IMPLEMENTED;
  return 0;
}

template <> int Tensor<int>::L2() const {
  NOT_IMPLEMENTED;
  return 0;
}

template <typename Dtype>
Dtype Tensor<Dtype>::L2() const {
  Dtype sumsq;
  const Dtype* data;
  if (!data_) { return 0; }
  switch (data_->head()) {
  case SyncedMemory::HEAD_AT_CPU:
    data = cpu_data();
    sumsq = cpu_dot(count_, data, data);
    break;
  case SyncedMemory::HEAD_AT_GPU:
  case SyncedMemory::SYNCED: //if data_ is already synced to cpu_data
#ifndef CPU_ONLY
    data = gpu_data();
    gpu_dot(count_, data, data, &sumsq);
#else
    NO_GPU;
#endif
    break;
  case SyncedMemory::UNINITIALIZED:
    return 0;
  default:
    std::cerr<< "Unknown SyncedMemory head state: " << data_->head()<<std::endl;
  }
  return sumsq;
}

/* This is a dangerous way filling tensor with values
 * Reshape() call must come before it to allocate mem */
template <typename Dtype>
void Tensor<Dtype>::reinit(const Dtype* source, unsigned int len) {
  CHECK_EQ(len,count_); 
  Copy(count_,source, static_cast<Dtype*>(data_->mutable_cpu_data()));
}


///* https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html */
//void Tensor<Dtype>::get(bp::array_t<double> xs) {
//    bp::buffer_info info = xs.request();
//    auto ptr = static_cast<double *>(info.ptr);
//    this->Reshape(info.shape);
//    this->reinit(info.ptr, info.size());
//} 


template <typename Dtype>
Tensor<Dtype> Tensor<Dtype>::operator[] (const unsigned int index) const {
  CHECK_LT(index, shape_[0]);
  Tensor<Dtype> result;
  result.Reshape(vector<unsigned int>(shape_.begin()+1,shape_.end()));
  result.reinit(cpu_data()+index*result.count(), result.count());
  return result;
}

template <typename Dtype>
Tensor<Dtype>& Tensor<Dtype>::operator= (const Tensor<Dtype>& otensor) {
  this->Reshape(otensor.shape());
  this->reinit(otensor.cpu_data(),otensor.count());
  return *this;
}

template <typename Dtype>
void Tensor<Dtype>::copy(const Tensor& source, unsigned int reshape) {
  if (source.count() != count_ || source.shape() != shape_) {
    CHECK (reshape);
    ReshapeLike(source);
  }

  switch (Axiom::mode()) {
  case Axiom::GPU:
    Copy(count_, source.gpu_data(), static_cast<Dtype*>(data_->mutable_gpu_data()));
    break;
  case Axiom::CPU:
    Copy(count_, source.cpu_data(), static_cast<Dtype*>(data_->mutable_cpu_data()));
    break;
  default:
    std::cerr<< "Unknown axiom mode."<<std::endl;
  }
}

INSTANTIATE_CLASS(Tensor);
template class Tensor<int>;
template class Tensor<unsigned int>;
} /* namespace axiom */

