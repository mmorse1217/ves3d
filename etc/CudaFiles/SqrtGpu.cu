#include <cuda.h>

#define BLOCK_HEIGHT 128

__global__
void sqrtKernel(const float *x_in, int length, float *x_out) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    x_out[idx] = sqrtf(x_in[idx]);
  }
}

void SqrtGpu(const float* x_in, int stride, int num_surfs, float *x_out) {
  int length = stride * num_surfs;
  int grid = length / BLOCK_HEIGHT + 1;
  sqrtKernel<<<grid, BLOCK_HEIGHT>>> (x_in, length, x_out);
}
