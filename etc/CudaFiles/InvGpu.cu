#include <cuda.h>

#define BLOCK_HEIGHT 128

__global__
void invKernel(const float *x_in, int length, float *x_out) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    x_out[idx] = 1.0F / x_in[idx];
  }
}

void InvGpu(const float* x_in, int stride, int num_surfs, float *x_out) {
  int length = stride * num_surfs;
  int grid = length / BLOCK_HEIGHT + 1;
  invKernel<<<grid, BLOCK_HEIGHT>>> (x_in, length, x_out);
}
