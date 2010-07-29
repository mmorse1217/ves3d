#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cuda.h>

extern void cuda_stokes(int, int, float*, float*, float*, float*);

void readfile(float *arr, char *fname, int size) {
  std::ifstream file(fname);
  if (file.is_open()) {
    int idx=0;
    while (idx < size) {
      file >> arr[idx++];
      fprintf(stdout, "%E\n", arr[idx - 1]);
    }
    fprintf(stdout, "Done\n");
    file.close();
  }
}


int main() {
  float *trg, *src, *den, *u, *q;
  int n = 1;
  int m = 544;
  int dim = 3;
  trg = (float*) malloc(sizeof(float) * dim * n);
  u = (float*) malloc(sizeof(float) * dim * n);
  src = (float*) malloc(sizeof(float) * dim * m * n);
  den = (float*) malloc(sizeof(float) * dim * m * n);
  readfile(trg, "targ_p12.txt", dim * n);
  readfile(src, "source_p12.txt", dim * m * n);
  readfile(den, "density_p12.txt", dim * m * n);

  cudaSetDevice(0);
  float *trg_dev, *src_dev, *den_dev, *u_dev, *q_dev;
  cudaMalloc((void**)&trg_dev, sizeof(float) * dim * n);
    cudaError_t C_E = cudaGetLastError ();
  fprintf (stderr, "%s\n", cudaGetErrorString (C_E));
  cudaMalloc((void**)&u_dev, sizeof(float) * dim * n);
  cudaMalloc((void**)&src_dev, sizeof(float) * dim * m * n);
  cudaMalloc((void**)&q_dev, sizeof(float) * dim * m * n);
  cudaMalloc((void**)&den_dev, sizeof(float) * dim * m * n);
  cudaMemcpy(trg_dev, trg, sizeof(float) * dim * n, cudaMemcpyHostToDevice);
    C_E = cudaGetLastError ();
  fprintf (stderr, "%s\n", cudaGetErrorString (C_E));
  cudaMemcpy(den_dev, den, sizeof(float) * dim * m * n, cudaMemcpyHostToDevice);
    C_E = cudaGetLastError ();
  fprintf (stderr, "%s\n", cudaGetErrorString (C_E));
  cudaMemcpy(src_dev, src, sizeof(float) * dim * m * n, cudaMemcpyHostToDevice);
    C_E = cudaGetLastError ();
  fprintf (stderr, "%s\n", cudaGetErrorString (C_E));

  cuda_stokes(m, n, trg_dev, src_dev, den_dev, u_dev);

  cudaMemcpy(u, u_dev, sizeof(float) * dim * n, cudaMemcpyDeviceToHost);

  for (int i = 0; i < dim * n; i++) {
    fprintf(stderr, "%E\n", u[i]);
  }
  return 0;
}
