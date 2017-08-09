#include "linalgops.hpp"

namespace axiom {

template <typename Dtype>
Dtype cpu_dot(const int n, const Dtype* x, const Dtype* y) {
  return cpu_strided_dot(n, x, 1, y, 1);
}
template
float cpu_dot<float>(const int n, const float* x, const float* y);
template
double cpu_dot<double>(const int n, const double* x, const double* y);

template <>
float cpu_strided_dot(const int n, const float* x, const int incx, const float* y, const int incy) {
  return cblas_sdot(n, x, incx, y, incy);
}

template <>
double cpu_strided_dot(const int n, const double* x, const int incx, const double* y, const int incy) {
  return cblas_ddot(n, x, incx, y, incy);
}

template <> //typename Dtype>
float cpu_asum<float>(const int n, const float* x) {
  return cblas_sasum(n, x, 1);
//  return 0;
}

template <> //typename Dtype>
double cpu_asum<double>(const int n, const double* x) {
  return cblas_dasum(n, x, 1);
//  return 0;
}
} /* namespace axiom */

