#include "linalgops.hpp"

namespace axiom {

template <>
void gpu_dot<float>(const int n, const float* x, const float* y,
    float* out) {
  CUBLAS_CHECK(cublasSdot(Axiom::cublas_handle(), n, x, 1, y, 1, out));
}
template <>
void gpu_dot<double>(const int n, const double* x, const double* y,
    double * out) {
  CUBLAS_CHECK(cublasDdot(Axiom::cublas_handle(), n, x, 1, y, 1, out));
}

template <>
void gpu_asum<float>(const int n, const float* x, float* y) {
  CUBLAS_CHECK(cublasSasum(Axiom::cublas_handle(), n, x, 1, y));
}
template <>
void gpu_asum<double>(const int n, const double* x, double* y) {
  CUBLAS_CHECK(cublasDasum(Axiom::cublas_handle(), n, x, 1, y));
}

} /* namespace axiom */
