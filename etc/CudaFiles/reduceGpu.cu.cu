#include <cuda.h>

#define BLOCK_HEIGHT 128

__global__
void reduceKernel(const float *x_in, const float *w_in, const float *q_in,
                   int stride, float *int_x_dw) {
  float threadSum = 0.0F;
  __shared__
  float sum[BLOCK_HEIGHT];

  int xOff = blockIdx.x * stride + threadIdx.x;
  int qOff = threadIdx.x;

  while(xOff < (blockIdx.x + 1)* stride) {
    threadSum = x_in[xOff] * w_in[xOff] * q_in[qOff];
    xOff += BLOCK_HEIGHT;
    qOff += BLOCK_HEIGHT;
  }
  sum[threadIdx.x] = threadSum;
  int redOff = 1;
  int redStride = 2;
  while(redOff != BLOCK_HEIGHT) {
    if (threadIdx.x % redStride == 0) {
      syncthreads();
      sum[threadIdx.x] += sum[threadIdx.x + redOff];
    }
    redOff = redStride;
    redStride *= 2;
  }
  if(threadIdx.x == 0) {
    int_x_dw[blockIdx.x] = sum[0];
  }
}


__global__
void reduceKernel(const float *w_in, const float *q_in,
                   int stride, float *int_x_dw) {
  float threadSum = 0.0F;
  __shared__
  float sum[BLOCK_HEIGHT];

  int xOff = blockIdx.x * stride + threadIdx.x;
  int qOff = threadIdx.x;

  while(xOff < (blockIdx.x + 1)* stride) {
    threadSum = w_in[xOff] * q_in[qOff];
    xOff += BLOCK_HEIGHT;
    qOff += BLOCK_HEIGHT;
  }
  sum[threadIdx.x] = threadSum;
  int redOff = 1;
  int redStride = 2;
  while(redOff != BLOCK_HEIGHT) {
    if (threadIdx.x % redStride == 0) {
      syncthreads();
      sum[threadIdx.x] += sum[threadIdx.x + redOff];
    }
    redOff = redStride;
    redStride *= 2;
  }
  if(threadIdx.x == 0) {
    int_x_dw[blockIdx.x] = sum[0];
  }
}


void ReduceGpu(const float *x_in, const float *w_in, const float *q_in,
                int stride, int num_surfs, float *int_x_dw) {
  int grid = num_surfs;
  if (x_in != NULL)
    reduceKernel<<<grid, BLOCK_HEIGHT>>> (x_in, w_in, q_in, stride, int_x_dw);
  else
    reduceKernel<<<grid, BLOCK_HEIGHT>>> (w_in, q_in, stride, int_x_dw);
}
