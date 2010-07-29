#include <cuda.h>
#include <stdio.h>

#define scalar float

#define BLOCK_HEIGHT 128


__global__
void shuffle(scalar *in, int m, int n, int dim, scalar *out) {
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


void cuda_shuffle(scalar *in, int m, int n, int dim, scalar *out) {
  int grid = n;
  shuffle<<<grid, BLOCK_HEIGHT>>> (in, m, n, dim, out);
}
