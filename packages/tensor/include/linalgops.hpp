#ifndef __LINALGOPS_H__
#define __LINALGOPS_H__

#include <stdint.h>
#include <math.h>
#include <cmath>
#ifdef __cplusplus
extern "C"
{
#endif
  #include <cblas.h>
#ifdef __cplusplus
}
#endif
#include "axiom.hpp"

namespace axiom {

template <typename Dtype>
Dtype cpu_strided_dot(const int n, const Dtype* x, const int incx,
                      const Dtype* y, const int incy);
template <typename Dtype>
Dtype cpu_dot(const int n, const Dtype* x, const Dtype* y);
/* Sum of absolute values of vector x */
template <typename Dtype>
Dtype cpu_asum(const int n, const Dtype* x);

#ifndef CPU_ONLY /* CPU + GPU */
  
template <typename Dtype>
void gpu_dot(const int n, const Dtype* x, const Dtype* y, Dtype* out);
/* Sum of absolute values of vector x */
template <typename Dtype>
void gpu_asum(const int n, const Dtype* x, Dtype* y);
  
#endif /* CPU + GPU */
}      /* namespace axiom */
#endif /* __LINALGOPS_H__ */
