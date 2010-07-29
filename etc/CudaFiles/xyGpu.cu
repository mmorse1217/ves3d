#include <cuda.h>

#define BLOCK_HEIGHT 128

__global__
void xyKernel(const float *x_in, const float *y_in, int length, float *xy_out) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    xy_out[idx] = x_in[idx] * y_in[idx];
  }
}

void xyGpu(const float* x_in, const float *y_in, int stride, int num_surfs, float *xy_out) {
  int length = stride * num_surfs;
  int grid = length / BLOCK_HEIGHT + 1;
  xyKernel<<<grid, BLOCK_HEIGHT>>> (x_in, y_in, length, xy_out);
}
