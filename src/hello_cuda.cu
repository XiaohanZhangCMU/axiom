#include "hello_cuda.hpp"

namespace axiom {

__global__ void saxpy(int n, float a, float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
//  printf("In GPU: n= %d: [%d,%d,%d]\n", n, blockIdx.x, blockDim.x, threadIdx.x);
}

int CudaAnimal::test_tensor_saxpy(void)
{
  Axiom::set_mode(Axiom::GPU); 

  unsigned int N = 1<<20; Tensor<float> x, y;
  std::vector<unsigned int> sp(1);
  sp[0] = N;
  x.Reshape(sp); y.Reshape(sp);
  float* xp = x.mutable_cpu_data();
  float* yp = y.mutable_cpu_data();
  for (int i = 0; i < N; i++) {
    xp[i] = 1.0f;
    yp[i] = 2.0f;
  }
  
  float* xg = x.mutable_gpu_data();
  float* yg = y.mutable_gpu_data();

  /* Perform SAXPY on 1M elements */
  saxpy<<<(N+255)/256, 256>>>(N, 2.0f, x.mutable_gpu_data(), y.mutable_gpu_data());

  const float* ypp = y.cpu_data();

  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
    maxError = max(maxError, abs(ypp[i]-4.0f));
  printf("Max error: %f\n", maxError);

  return 0;
}
int CudaAnimal::test_saxpy(void)
{
  int N = 1<<20;
  float *x, *y, *d_x, *d_y;
  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));

  cudaMalloc(&d_x, N*sizeof(float));
  cudaMalloc(&d_y, N*sizeof(float));

  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }
//  std::cout<<"This is a test"<<std::endl;
  cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice);

  // Perform SAXPY on 1M elements
  saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y);

  cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);

  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
    maxError = max(maxError, abs(y[i]-4.0f));
  printf("Max error: %f\n", maxError);

  cudaFree(d_x);
  cudaFree(d_y);
  free(x);
  free(y);

  return 0;
}

}
