#include "Device.h"

extern const Device<CPU> device(0);

int main(int argc, char **argv)
{
  int k = 4096;
  int n = 4096;

  size_t sa = n*n;
  size_t sb = n*k;
  size_t sc = n*k;

  float *A_h = (float*) malloc(sa * sizeof(float));
  float *B_h = (float*) malloc(sb * sizeof(float));

  float *A = (float*) device.Malloc(sa * sizeof(float));
  float *B = (float*) device.Malloc(sb * sizeof(float));
  float *C = (float*) device.Malloc(sc * sizeof(float));

  float alpha = 1.0;
  float beta = 0.0;

  for(int ii=0;ii<n;++ii)
    for(int jj=0;jj<n;++jj)
      A_h[ii*n+jj] = drand48();

  for(int ii=0;ii<n;++ii)
    for(int jj=0;jj<k;++jj)
      B_h[ii*k+jj] = drand48();

  device.Memcpy(A,A_h,sa,MemcpyHostToDevice);
  device.Memcpy(B,B_h,sa,MemcpyHostToDevice);

  int i_max = 1;
  for(int ii=0;ii<i_max;++ii)
  {
      device.gemm("N", "N", &n, &k, &n, &alpha, A, &n, B, &n, &beta, C, &n);
  }
  
  PROFILEREPORT(SortFlopRate);

  device.Memcpy(A_h,C,sa,MemcpyDeviceToHost);

  device.Free(A);
  device.Free(B);
  device.Free(C);

  free(A_h);
  free(B_h);


}
