#include <cuda.h>

#define BLOCK_HEIGHT 128

__global__
void avpwKernel(const float *a_in, const float *v_in, const float *w_in,
                 int size, int length, float *avpw_out) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    avpw_out[idx] = a_in[idx / (size)] * v_in[idx] + w_in[idx];
  }
}


void avpwGpu(const float *a_in, const float *v_in, const float *w_in,
              int stride, int num_surfs, float *avpw_out) {
  dim3 grid;
  grid.x = num_surfs * stride * 3 / BLOCK_HEIGHT + 1;
  avpwKernel<<<grid, BLOCK_HEIGHT>>> (a_in, v_in, w_in, stride * 3, num_surfs * stride * 3, avpw_out);
}
