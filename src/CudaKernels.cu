#include "CudaKernels.h"
#include <stdio.h>
#include <cuda.h>

#define DIM 3

#define BLOCK_HEIGHT 128

#define DOT_FLOPS 5
#define CROSS_FLOPS 9
#define AXPY_FLOPS 6
#define XDY_FLOPS 3

void CE() {
  cudaError_t ce = cudaGetLastError();
  printf("%s\n", cudaGetErrorString(ce));
}

__global__
void dotProdKernel(const float* a, const float* b, int stride, int num_surfs, float* aDb) {
  unsigned int blockOff, resOff, off, res;
  float aXReg, aYReg, aZReg, bXReg, bYReg, bZReg, dotProd;
  resOff = blockIdx.x * stride;
  blockOff = resOff * DIM;

  int numChunkLoop = stride / BLOCK_HEIGHT;

  off = blockOff + threadIdx.x;
  res = resOff + threadIdx.x;
  
  for (int chunk = 0; chunk < numChunkLoop; chunk++) {
    aXReg = a[off];
    aYReg = a[off + stride];
    aZReg = a[off + stride + stride];
    bXReg = b[off];
    bYReg = b[off + stride];
    bZReg = b[off + stride + stride];
    dotProd = aXReg * bXReg + aYReg * bYReg + aZReg * bZReg;
    aDb[res] = dotProd;
    off += BLOCK_HEIGHT;
    res += BLOCK_HEIGHT;
  }

  if (off < blockOff + stride) {
    aXReg = a[off];
    aYReg = a[off + stride];
    aZReg = a[off + stride + stride];
    bXReg = b[off];
    bYReg = b[off + stride];
    bZReg = b[off + stride + stride];
    dotProd = aXReg * bXReg + aYReg * bYReg + aZReg * bZReg;
    aDb[res] = dotProd;
  }
}

void DotProductGpu(const float *a, const float *b, int stride, int num_surfs, float *aDb)
{
#ifdef GPU_PROF
    float kernelTime, flops, gflopRate;
    flops = DOT_FLOPS * stride * num_surfs;
    cudaEvent_t kernelStart, kernelStop;
    cudaEventCreate(&kernelStart);
    cudaEventCreate(&kernelStop);
    cudaEventRecord(kernelStart, 0);
#endif

    int GridDim = num_surfs;

  dotProdKernel<<<GridDim, BLOCK_HEIGHT>>>
      (a, b, stride, num_surfs, aDb);

  cudaThreadSynchronize();
#ifdef GPU_PROF
  cudaEventRecord(kernelStop, 0);
  cudaEventSynchronize(kernelStop);
  cudaEventElapsedTime(&kernelTime, kernelStart, kernelStop);
  gflopRate = (flops / 1e9) / (kernelTime / 1e3);
  fprintf(stderr, "Dot product kernel took %f ms @ %f Gflops\n", kernelTime, gflopRate);
#endif
}

__global__
void crossProdKernel(const float* a, const float* b, int stride, int num_surfs, float* aCb) {
  unsigned int blockOff, off;
  float aXReg, aYReg, aZReg, bXReg, bYReg, bZReg, aCbXReg, aCbYReg, aCbZReg;
  blockOff = blockIdx.x * stride * DIM;

  int numChunkLoop = stride / BLOCK_HEIGHT;

  off = blockOff + threadIdx.x;
  
  for (int chunk = 0; chunk < numChunkLoop; chunk++) {
    aXReg = a[off];
    aYReg = a[off + stride];
    aZReg = a[off + stride + stride];
    bXReg = b[off];
    bYReg = b[off + stride];
    bZReg = b[off + stride + stride];

    aCbXReg = aYReg * bZReg - aZReg * bYReg;
    aCbYReg = aZReg * bXReg - aXReg * bZReg;
    aCbZReg = aXReg * bYReg - aYReg * bXReg;

    aCb[off] = aCbXReg;
    aCb[off + stride] = aCbYReg;
    aCb[off + stride + stride] = aCbZReg;

    off += BLOCK_HEIGHT;
  }

  if (off < blockOff + stride) {
    aXReg = a[off];
    aYReg = a[off + stride];
    aZReg = a[off + stride + stride];
    bXReg = b[off];
    bYReg = b[off + stride];
    bZReg = b[off + stride + stride];

    aCbXReg = aYReg * bZReg - aZReg * bYReg;
    aCbYReg = aZReg * bXReg - aXReg * bZReg;
    aCbZReg = aXReg * bYReg - aYReg * bXReg;

    aCb[off] = aCbXReg;
    aCb[off + stride] = aCbYReg;
    aCb[off + stride + stride] = aCbZReg;
  }
}

void CrossProductGpu(const float *a, const float *b, int stride, int num_surfs, float *aCb) 
{
#ifdef GPU_PROF
  float kernelTime, flops, gflopRate;
  flops = CROSS_FLOPS * stride * num_surfs;
  cudaEvent_t kernelStart, kernelStop;
  cudaEventCreate(&kernelStart);
  cudaEventCreate(&kernelStop);
  cudaEventRecord(kernelStart, 0);
#endif

  int GridDim = num_surfs;
  crossProdKernel<<<GridDim, BLOCK_HEIGHT>>>
      (a, b, stride, num_surfs, aCb);

  cudaThreadSynchronize();
#ifdef GPU_PROF
  cudaEventRecord(kernelStop, 0);
  cudaEventSynchronize(kernelStop);
  cudaEventElapsedTime(&kernelTime, kernelStart, kernelStop);
  gflopRate = (flops / 1e9) / (kernelTime / 1e3);
  fprintf(stderr, "Cross product kernel took %f ms @ %f Gflops\n", kernelTime, gflopRate);
#endif
}

__global__
void xvpwKernel(const float *x, const float *a, const float *y, int stride, int num_surfs, float *AxPy) {
  unsigned int blockOff, scalOff, off, scal;
  float xXReg, xYReg, xZReg, aReg, yXReg, yYReg, yZReg, AxPyXReg, AxPyYReg, AxPyZReg;
  scalOff = blockIdx.x * stride;
  blockOff = scalOff * DIM;

  int numChunkLoop = stride / BLOCK_HEIGHT;

  off = blockOff + threadIdx.x;
  scal = scalOff + threadIdx.x;
  
  for (int chunk = 0; chunk < numChunkLoop; chunk++) {
    xXReg = x[off];
    xYReg = x[off + stride];
    xZReg = x[off + stride + stride];
    yXReg = y[off];
    yYReg = y[off + stride];
    yZReg = y[off + stride + stride];
    aReg = a[scal];
    

    AxPyXReg = aReg * xXReg + yXReg;
    AxPyYReg = aReg * xYReg + yYReg;
    AxPyZReg = aReg * xZReg + yZReg;

    AxPy[off] = AxPyXReg;
    AxPy[off + stride] = AxPyYReg;
    AxPy[off + stride + stride] = AxPyZReg;

    off += BLOCK_HEIGHT;
    scal += BLOCK_HEIGHT;
  }

  if (off < blockOff + stride) {
    xXReg = x[off];
    xYReg = x[off + stride];
    xZReg = x[off + stride + stride];
    yXReg = y[off];
    yYReg = y[off + stride];
    yZReg = y[off + stride + stride];
    aReg = a[scal];
    

    AxPyXReg = aReg * xXReg + yXReg;
    AxPyYReg = aReg * xYReg + yYReg;
    AxPyZReg = aReg * xZReg + yZReg;

    AxPy[off] = AxPyXReg;
    AxPy[off + stride] = AxPyYReg;
    AxPy[off + stride + stride] = AxPyZReg;
  }

}

void xvpwGpu(const float *x, const float *a, const float *y, int stride, int num_surfs, float *AxPy) {

#ifdef GPU_PROF
  float kernelTime, flops, gflopRate;
  flops = AXPY_FLOPS * stride * num_surfs;
  cudaEvent_t kernelStart, kernelStop;
  cudaEventCreate(&kernelStart);
  cudaEventCreate(&kernelStop);
  cudaEventRecord(kernelStart, 0);
#endif

  int GridDim = num_surfs;
  xvpwKernel<<<GridDim, BLOCK_HEIGHT>>>
      (x, a, y, stride, num_surfs, AxPy);

  cudaThreadSynchronize();

#ifdef GPU_PROF
  cudaEventRecord(kernelStop, 0);
  cudaEventSynchronize(kernelStop);
  cudaEventElapsedTime(&kernelTime, kernelStart, kernelStop);
  gflopRate = (flops / 1e9) / (kernelTime / 1e3);
  fprintf(stderr, "AxPy kernel took %f ms @ %f Gflops\n", kernelTime, gflopRate);
#endif
}

__global__
void xvpbKernel(const float *x, const float *a, float y, int stride, int num_surfs, float *AxPy) {
  unsigned int blockOff, scalOff, off, scal;
  float xXReg, xYReg, xZReg, aReg, AxPyXReg, AxPyYReg, AxPyZReg;
  scalOff = blockIdx.x * stride;
  blockOff = scalOff * DIM;

  int numChunkLoop = stride / BLOCK_HEIGHT;

  off = blockOff + threadIdx.x;
  scal = scalOff + threadIdx.x;
  
  for (int chunk = 0; chunk < numChunkLoop; chunk++) {
    xXReg = x[off];
    xYReg = x[off + stride];
    xZReg = x[off + stride + stride];
    aReg = a[scal];
    

    AxPyXReg = aReg * xXReg + y;
    AxPyYReg = aReg * xYReg + y;
    AxPyZReg = aReg * xZReg + y;

    AxPy[off] = AxPyXReg;
    AxPy[off + stride] = AxPyYReg;
    AxPy[off + stride + stride] = AxPyZReg;

    off += BLOCK_HEIGHT;
    scal += BLOCK_HEIGHT;
  }

  if (off < blockOff + stride) {
    xXReg = x[off];
    xYReg = x[off + stride];
    xZReg = x[off + stride + stride];
    aReg = a[scal];
    

    AxPyXReg = aReg * xXReg + y;
    AxPyYReg = aReg * xYReg + y;
    AxPyZReg = aReg * xZReg + y;

    AxPy[off] = AxPyXReg;
    AxPy[off + stride] = AxPyYReg;
    AxPy[off + stride + stride] = AxPyZReg;
  }
}

void xvpbGpu(const float *x, const float *a, float y, int stride, int num_surfs, float *AxPy) {

#ifdef GPU_PROF
  float kernelTime, flops, gflopRate;
  flops = AXPY_FLOPS * stride * num_surfs;
  cudaEvent_t kernelStart, kernelStop;
  cudaEventCreate(&kernelStart);
  cudaEventCreate(&kernelStop);
  cudaEventRecord(kernelStart, 0);
#endif

  int GridDim = num_surfs;
  xvpbKernel<<<GridDim, BLOCK_HEIGHT>>>
      (x, a, y, stride, num_surfs, AxPy);

  cudaThreadSynchronize();

#ifdef GPU_PROF
  cudaEventRecord(kernelStop, 0);
  cudaEventSynchronize(kernelStop);
  cudaEventElapsedTime(&kernelTime, kernelStart, kernelStop);
  gflopRate = (flops / 1e9) / (kernelTime / 1e3);
  fprintf(stderr, "AxPy kernel took %f ms @ %f Gflops\n", kernelTime, gflopRate);
#endif
}

__global__
void uyInvKernel(const float *x, const float *y, int stride, int num_surfs, float *xDy) {
  unsigned int blockOff, scalOff, off, scal;
  float xXReg, xYReg, xZReg, yReg, xDyXReg, xDyYReg, xDyZReg;
  scalOff = blockIdx.x * stride;
  blockOff = scalOff * DIM;

  int numChunkLoop = stride / BLOCK_HEIGHT;

  off = blockOff + threadIdx.x;
  scal = scalOff + threadIdx.x;
  
  for (int chunk = 0; chunk < numChunkLoop; chunk++) {
    xXReg = x[off];
    xYReg = x[off + stride];
    xZReg = x[off + stride + stride];
    yReg = y[scal];
    

    xDyXReg = xXReg / yReg;
    xDyYReg = xYReg / yReg;
    xDyZReg = xZReg / yReg;

    xDy[off] = xDyXReg;
    xDy[off + stride] = xDyYReg;
    xDy[off + stride + stride] = xDyZReg;

    off += BLOCK_HEIGHT;
    scal += BLOCK_HEIGHT;
  }

  if (off < blockOff + stride) {
    xXReg = x[off];
    xYReg = x[off + stride];
    xZReg = x[off + stride + stride];
    yReg = y[scal];
    

    xDyXReg = xXReg / yReg;
    xDyYReg = xYReg / yReg;
    xDyZReg = xZReg / yReg;

    xDy[off] = xDyXReg;
    xDy[off + stride] = xDyYReg;
    xDy[off + stride + stride] = xDyZReg;
  }
}

void uyInvGpu(const float *x, const float *y, int stride, int num_surfs, float *xDy) {

#ifdef GPU_PROF
  float kernelTime, flops, gflopRate;
  flops = XDY_FLOPS * stride * num_surfs;
  cudaEvent_t kernelStart, kernelStop;
  cudaEventCreate(&kernelStart);
  cudaEventCreate(&kernelStop);
  cudaEventRecord(kernelStart, 0);
#endif

  int GridDim = num_surfs;
  uyInvKernel<<<GridDim, BLOCK_HEIGHT>>>
      (x, y, stride, num_surfs, xDy);

  cudaThreadSynchronize();

#ifdef GPU_PROF
  cudaEventRecord(kernelStop, 0);
  cudaEventSynchronize(kernelStop);
  cudaEventElapsedTime(&kernelTime, kernelStart, kernelStop);
  gflopRate = (flops / 1e9) / (kernelTime / 1e3);
  fprintf(stderr, "xDy kernel took %f ms @ %f Gflops\n", kernelTime, gflopRate);
#endif
}

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
  cudaThreadSynchronize();
}

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
  cudaThreadSynchronize();
}

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
  cudaThreadSynchronize();
}

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
  cudaThreadSynchronize();
}

__global__
void reduceKernel(const float *x_in, const float *w_in, const float *q_in,
                   int stride, float *int_x_dw) {
  float threadSum = 0.0F;
  __shared__
  float sum[BLOCK_HEIGHT];

  int xOff = blockIdx.x * stride + threadIdx.x;
  int qOff = threadIdx.x;

  while(xOff < (blockIdx.x + 1)* stride) {
    threadSum += x_in[xOff] * w_in[xOff] * q_in[qOff];
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
    threadSum += w_in[xOff] * q_in[qOff];
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
  cudaThreadSynchronize();
}

void CircShiftGpu(const float *arr_in, int n_vecs, int vec_length, int shift, float *arr_out) {
  shift = shift % vec_length;
  if (shift < 0) {
    shift += vec_length;
  }
  int base_in, base_out;
  for (int ii = 0; ii < n_vecs; ii++) {
    base_out = ii * vec_length;
    base_in = base_out + vec_length - shift;
    cudaMemcpy(arr_out + base_out, arr_in + base_in, sizeof(float) * shift, cudaMemcpyDeviceToDevice);
    base_in = base_out;
    base_out += shift;
    cudaMemcpy(arr_out + base_out, arr_in + base_in, sizeof(float) * (vec_length - shift), cudaMemcpyDeviceToDevice);
  }
  cudaThreadSynchronize();
}

__global__
void axpyKernel(float a, const float *x, const float *y, int length, float *axpy) {

  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    axpy[idx] = a * x[idx] + y[idx];
  }
}

void axpyGpu(float a, const float* x_in, const float *y_in, int stride, int num_surfs, float *axpy_out) {
  int length = stride * num_surfs;
  int grid = length / BLOCK_HEIGHT + 1;
  axpyKernel<<<grid, BLOCK_HEIGHT>>> (a, x_in, y_in, length, axpy_out);
}

__global__
void axpbKernel(float a, const float *x, const float b, int length, float *axpy) {

  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    axpy[idx] = a * x[idx] + b;
  }
}

void axpbGpu(float a, const float* x_in, float b, int stride, int num_surfs, float *axpb_out) {
  int length = stride * num_surfs;
  int grid = length / BLOCK_HEIGHT + 1;
  axpbKernel<<<grid, BLOCK_HEIGHT>>> (a, x_in, b, length, axpb_out);
  cudaThreadSynchronize();
}

__global__
void shuffle(float *in, int m, int n, int dim, float *out) {
  int sub, add;
  int block_off = blockIdx.x * dim * m;

  in += block_off;
  out += block_off;

  int thread_off = threadIdx.x;
  int out_off;

  while(thread_off < dim * m) {
    int f = thread_off / m;
    sub = f * m;
    add = f;
    out_off = (thread_off - sub) * dim + add;
    out[out_off] = in[thread_off];

    thread_off += BLOCK_HEIGHT;
  }
}


void cuda_shuffle(float *in, int m, int n, int dim, float *out) {
  int grid = n;
  shuffle<<<grid, BLOCK_HEIGHT>>> (in, m, n, dim, out);
  cudaThreadSynchronize();
}

///
#define PI_8I 0.0397887358F

__global__
void stokes(int m, int n, int t_head, const float *T, const float *S, const float *D, float *U, const float *Q) {
  float3 trg_reg;
  float3 src_reg;
  float3 pot_reg;
  float3 dis_reg;
  float3 u_reg;
  __shared__
  float3 u_sh[BLOCK_HEIGHT];

  int t_off = blockIdx.x * 3 * m + t_head + blockIdx.y;

  trg_reg.x = T[t_off];
  trg_reg.y = T[m + t_off];
  trg_reg.z = T[m + m + t_off];

  u_reg = make_float3(0.0, 0.0, 0.0);

  int block_off = blockIdx.x * 3 * m + threadIdx.x;
  int s = threadIdx.x;

  while(block_off < blockIdx.x * 3 * m + m) {
    src_reg.x = S[block_off];
    src_reg.y = S[block_off + m];
    src_reg.z = S[block_off + m + m];

    pot_reg.x = D[block_off] * Q[s];
    pot_reg.y = D[block_off + m] * Q[s];
    pot_reg.z = D[block_off + m + m] * Q[s];

    dis_reg.x = src_reg.x - trg_reg.x;
    dis_reg.y = src_reg.y - trg_reg.y;
    dis_reg.z = src_reg.z - trg_reg.z;

    float inv_r = rsqrtf(dis_reg.x * dis_reg.x + dis_reg.y * dis_reg.y
                          + dis_reg.z * dis_reg.z);

    inv_r = inv_r + (inv_r-inv_r);
    inv_r = fmaxf(inv_r,0.0F);
    
    float tmp_scal = (dis_reg.x * pot_reg.x + dis_reg.y * pot_reg.y
                       + dis_reg.z * pot_reg.z) * inv_r * inv_r;
    pot_reg.x += tmp_scal * dis_reg.x;
    pot_reg.y += tmp_scal * dis_reg.y;
    pot_reg.z += tmp_scal * dis_reg.z;

    u_reg.x += pot_reg.x * inv_r;
    u_reg.y += pot_reg.y * inv_r;
    u_reg.z += pot_reg.z * inv_r;

    block_off += BLOCK_HEIGHT;
    s += BLOCK_HEIGHT;
  }

  u_sh[threadIdx.x].x = u_reg.x;
  u_sh[threadIdx.x].y = u_reg.y;
  u_sh[threadIdx.x].z = u_reg.z;

  int off = 1;
  int stride = 2;
  while (off != BLOCK_HEIGHT) {
    if (threadIdx.x % stride == 0) {
      syncthreads();
      u_sh[threadIdx.x].x += u_sh[threadIdx.x + off].x;
      u_sh[threadIdx.x].y += u_sh[threadIdx.x + off].y;
      u_sh[threadIdx.x].z += u_sh[threadIdx.x + off].z;
    }
    off = stride;
    stride *= 2;
  }
  if (threadIdx.x == 0) {
    U[t_off] = u_sh[0].x * PI_8I;
    U[m + t_off] = u_sh[0].y * PI_8I;
    U[m + m + t_off] = u_sh[0].z * PI_8I;
  }

}

__global__
void stokes(int m, int n, int t_head, const float *T, const float *S, const float *D, float *U) {
  float3 trg_reg;
  float3 src_reg;
  float3 pot_reg;
  float3 dis_reg;
  float3 u_reg;
  __shared__
  float3 u_sh[BLOCK_HEIGHT];

  int t_off = blockIdx.x * 3 * m + t_head + blockIdx.y;

  trg_reg.x = T[t_off];
  trg_reg.y = T[m + t_off];
  trg_reg.z = T[m + m + t_off];

  u_reg = make_float3(0.0, 0.0, 0.0);

  int block_off = blockIdx.x * 3 * m + threadIdx.x;
  int s = threadIdx.x;

  while(block_off < blockIdx.x * 3 * m + m) {
    src_reg.x = S[block_off];
    src_reg.y = S[block_off + m];
    src_reg.z = S[block_off + m + m];

    pot_reg.x = D[block_off];
    pot_reg.y = D[block_off + m];
    pot_reg.z = D[block_off + m + m];

    dis_reg.x = src_reg.x - trg_reg.x;
    dis_reg.y = src_reg.y - trg_reg.y;
    dis_reg.z = src_reg.z - trg_reg.z;

    float inv_r = rsqrtf(dis_reg.x * dis_reg.x + dis_reg.y * dis_reg.y
                          + dis_reg.z * dis_reg.z);

    inv_r = inv_r + (inv_r-inv_r);
    inv_r = fmaxf(inv_r,0.0F);
    
    float tmp_scal = (dis_reg.x * pot_reg.x + dis_reg.y * pot_reg.y
                       + dis_reg.z * pot_reg.z) * inv_r * inv_r;
    pot_reg.x += tmp_scal * dis_reg.x;
    pot_reg.y += tmp_scal * dis_reg.y;
    pot_reg.z += tmp_scal * dis_reg.z;

    u_reg.x += pot_reg.x * inv_r;
    u_reg.y += pot_reg.y * inv_r;
    u_reg.z += pot_reg.z * inv_r;

    block_off += BLOCK_HEIGHT;
    s += BLOCK_HEIGHT;
  }

  u_sh[threadIdx.x].x = u_reg.x;
  u_sh[threadIdx.x].y = u_reg.y;
  u_sh[threadIdx.x].z = u_reg.z;

  int off = 1;
  int stride = 2;
  while (off != BLOCK_HEIGHT) {
    if (threadIdx.x % stride == 0) {
      syncthreads();
      u_sh[threadIdx.x].x += u_sh[threadIdx.x + off].x;
      u_sh[threadIdx.x].y += u_sh[threadIdx.x + off].y;
      u_sh[threadIdx.x].z += u_sh[threadIdx.x + off].z;
    }
    off = stride;
    stride *= 2;
  }
  if (threadIdx.x == 0) {
    U[t_off] = u_sh[0].x * PI_8I;
    U[m + t_off] = u_sh[0].y * PI_8I;
    U[m + m + t_off] = u_sh[0].z * PI_8I;
  }

}


void cuda_stokes(int m, int n, int t_head, int t_tail, const float *T, const float *S, const float *D, float *U, const float *Q) {
  dim3 grid;
  grid.x = n;
  grid.y = t_tail - t_head;

  if (Q != NULL)
      stokes<<<grid, BLOCK_HEIGHT>>> (m, n, t_head, T, S, D, U, Q);
  else
      stokes<<<grid, BLOCK_HEIGHT>>> (m, n, t_head, T, S, D, U);
  cudaThreadSynchronize();
}

void ResampleGpu(int p, int n_funs, int q, const float *shc_p, float *shc_q) {

  float *out_deb = shc_q;
  int leg_order = p + 1;
  int new_leg_order = q + 1;
  int min_leg_order = (leg_order < new_leg_order) ? leg_order : new_leg_order;

  for(int v = 0; v < n_funs; v++) {
    cudaMemcpy(shc_q, shc_p, sizeof(float) * min_leg_order, cudaMemcpyDeviceToDevice);
    shc_q += min_leg_order;
    shc_p += min_leg_order;
    if (new_leg_order > leg_order) {
      cudaMemset(shc_q, 0, sizeof(float) * (new_leg_order - leg_order));
      shc_q += (new_leg_order - leg_order);
    }
    if (leg_order > new_leg_order)
      shc_p += (leg_order - new_leg_order);
  }
  leg_order--;
  new_leg_order--;
  min_leg_order--;

  for(; min_leg_order > 1; min_leg_order--, leg_order--, new_leg_order--) {
    for(int v = 0; v < n_funs; v++) {
      cudaMemcpy(shc_q, shc_p, sizeof(float) * min_leg_order, cudaMemcpyDeviceToDevice);
      shc_q += min_leg_order;
      shc_p += min_leg_order;
      if (new_leg_order > leg_order) {
        cudaMemset(shc_q, 0, sizeof(float) * (new_leg_order - leg_order));
        shc_q += (new_leg_order - leg_order);
      }
      if (leg_order > new_leg_order)
        shc_p += (leg_order - new_leg_order);
    }
    for(int v = 0; v < n_funs; v++) {
      cudaMemcpy(shc_q, shc_p, sizeof(float) * min_leg_order, cudaMemcpyDeviceToDevice);
      shc_q += min_leg_order;
      shc_p += min_leg_order;
      if (new_leg_order > leg_order) {
        cudaMemset(shc_q, 0, sizeof(float) * (new_leg_order - leg_order));
        shc_q += (new_leg_order - leg_order);
      }
      if (leg_order > new_leg_order)
        shc_p += (leg_order - new_leg_order);
    }
  }

  for(int v = 0; v < n_funs; v++) {
    cudaMemcpy(shc_q, shc_p, sizeof(float) * min_leg_order, cudaMemcpyDeviceToDevice);
    shc_q += min_leg_order;
    shc_p += min_leg_order;
    if (new_leg_order > leg_order) {
      cudaMemset(shc_q, 0, sizeof(float) * (new_leg_order - leg_order));
      shc_q += (new_leg_order - leg_order);
    }
    if (leg_order > new_leg_order)
      shc_p += (leg_order - new_leg_order);
  }

  leg_order--;
  new_leg_order--;
  min_leg_order--;

  float *outputs_end = out_deb + n_funs * q * (q + 2);
  if (shc_q < outputs_end) {
    cudaMemset(shc_q, 0, sizeof(float) * (outputs_end - shc_q));
  }
  cudaThreadSynchronize();
}

__global__
void xyMKernel(const float *x_in, const float *y_in, int length, int stride, float *xy_out) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    xy_out[idx] = x_in[idx] * y_in[idx % stride];
  }
}

void xyMGpu(const float* x_in, const float *y_in, int stride, int num_surfs, float *xy_out) {
  int length = stride * num_surfs;
  int grid = length / BLOCK_HEIGHT + 1;
  xyMKernel<<<grid, BLOCK_HEIGHT>>> (x_in, y_in, length, stride, xy_out);
  cudaThreadSynchronize();
}

void ScaleFreqsGpu(int p, int n_funs, const float *shc_in, const float *alpha, float *shc_out) {
    int leg_order = p+1;
    xyMGpu(shc_in, alpha, leg_order, n_funs, shc_out);
    alpha += leg_order;
    shc_in += n_funs * leg_order;
    shc_out += n_funs * leg_order;
    leg_order--;

    // process remaining frequencies except the last cosine
    for (; leg_order>1; leg_order--) 
    {
        // first process cosine
        xyMGpu(shc_in, alpha, leg_order, n_funs, shc_out);
        alpha += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
        
        // then process sine
        xyMGpu(shc_in, alpha, leg_order, n_funs, shc_out);
        alpha += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
    }
    
    // process last cosine
    xyMGpu(shc_in, alpha, leg_order, n_funs, shc_out);
    cudaThreadSynchronize();
}

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
  cudaThreadSynchronize();
}


__global__
void reduceMaxKernel(float *in, int n) {
  __shared__ float sdata[BLOCK_HEIGHT];
  int idx = blockIdx.x * BLOCK_HEIGHT + threadIdx.x;
  if (idx < n)
    sdata[threadIdx.x] = in[idx];
  else
    sdata[threadIdx.x] = -1e9;

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
    reduceMaxKernel<<<grid, BLOCK_HEIGHT>>> (in, n);
    n /= BLOCK_HEIGHT;
  }
  float max;
  cudaMemcpy(&max, in, sizeof(float), cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
  return max;
}
