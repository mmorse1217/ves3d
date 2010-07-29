#include <cuda.h>
#include <stdio.h>

extern void cuda_shuffle(float *in, int m, int n);

int main() {
  float *in;
  int m = 24;
  int n = 100;

  in = (float*) malloc(sizeof(float) * 3 * m * n);

  for (int block = 0; block < n; block++) {
    for (int dim = 0; dim < 3; dim++) {
      for (int idx = 0; idx < m; idx++) {
        in[block * 3 * m + dim * m + idx] = block * 1000 + dim * 100 + idx;
      }
    }
  }
  float *in_dev;
  cudaMalloc(&in_dev, sizeof(float) * 3 * m * n);

  cudaMemcpy(in_dev, in, sizeof(float) * 3 * m * n, cudaMemcpyHostToDevice);

  cuda_shuffle(in_dev, m, n);

  cudaMemcpy(in, in_dev, sizeof(float) * 3 * m * n, cudaMemcpyDeviceToHost);

  for (int i = 0; i < 3 * m * n; i++) {
    fprintf(stdout, "%4f\n", in[i]);
  }
  return 0;
}
