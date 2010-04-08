#include "DeviceCPU.h"

int main(int argc, char *arvg[])
{
  int n = 1000;
  int k = 4000;
  
  size_t sa = n*n;
  size_t sb = n*k;
  size_t sc = n*k;

  DeviceCPU<float> cpu;

  float *A = cpu.Malloc(sa);
  float *B = cpu.Malloc(sb);
  float *C = cpu.Malloc(sc);

  float alpha = 1.0;
  float beta = 0.0;

  for(int ii=0;ii<n;++ii)
    for(int jj=0;jj<n;++jj)
      A[ii*n+jj] = drand48();

  for(int ii=0;ii<n;++ii)
    for(int jj=0;jj<k;++jj)
      B[ii*k+jj] = drand48();

  double ss = get_seconds();
  for(int ii=0;ii<100;++ii)
    cpu.gemm("N", "N", &n, &k, &n, &alpha, A, &n, B, &n, &beta, C, &n);

  ss = get_seconds()-ss;
  cout<<"Time :"<<ss<<endl;

  cpu.Free(A);
  cpu.Free(B);
  cpu.Free(C);
}
