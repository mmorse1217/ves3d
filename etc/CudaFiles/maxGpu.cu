#include <cuda.h>

#define BLOCK_HEIGHT 128

__global__
void reduceKernel(float *in, int n) {
  __shared__ float sdata[BLOCK_HEIGHT];
  int idx = blockIdx.x * BLOCK_HEIGHT + threadIdx.x;
  if (idx < n)
    sdata[threadIdx.x] = in[idx];
  else
    sdata[threadIdx.x] = 1e-9;

  int redOff = 1;
  int redStride = 2;
  while(redOff != BLOCK_HEIGHT) {
    if (threadIdx.x % redStride == 0) {
      syncthreads();
      sdata[threadIdx.x] = fmaxf(sdata[threadIdx.x], sdata[threadIdx.x + redOff]);
    }
    redOff = redStride;
    redStride *= 2;
  }
  if(threadIdx.x == 0) {
    in[blockIdx.x] = sdata[0];
  }
}

float maxGpu(float *in, int n) {
  while(n > 0) {
    int grid = n / BLOCK_HEIGHT + 1;
    reduceKernel<<<grid, BLOCK_HEIGHT>>> (in, n);
    n /= BLOCK_HEIGHT;
  }
  float max;
  cudaMemcpy(&max, in, sizeof(float), cudaMemcpyDeviceToHost);
  return max;
}
