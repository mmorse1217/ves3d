#include <cuda.h>

#define BLOCK_HEIGHT 128

__global__
void xyInvKernel(const float *x_in, const float *y_in, int length, float *xDy_out) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    xDy_out[idx] = x_in[idx] / y_in[idx];
  }
}

void xyInvGpu(const float* x_in, const float *y_in, int stride, int num_surfs, float *xDy_out) {
  int length = stride * num_surfs;
  int grid = length / BLOCK_HEIGHT + 1;
  xyInvKernel<<<grid, BLOCK_HEIGHT>>> (x_in, y_in, length, xDy_out);
}
