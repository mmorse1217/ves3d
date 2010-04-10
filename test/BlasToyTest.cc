//#include "DeviceCPU.h"
#include "DeviceGPU.h"

int main(int argc, char *arvg[])
{
  int n = 4096*2;
  int k = 4096*2;
  
  size_t sa = n*n;
  size_t sb = n*k;
  size_t sc = n*k;

  DeviceGPU<float> device;

  float *A_h = (float*) malloc(sa * sizeof(float));
  float *B_h = (float*) malloc(sb * sizeof(float));

  float *A = device.Malloc(sa);
  float *B = device.Malloc(sb);
  float *C = device.Malloc(sc);

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
  double ss = get_seconds();
  for(int ii=0;ii<i_max;++ii)
  {
      device.gemm("N", "N", &n, &k, &n, &alpha, A, &n, B, &n, &beta, C, &n);
      cudaThreadSynchronize();
  }
  ss = get_seconds()-ss;
  cout<<"Time : "<<ss<<endl;
  cout<<"GFLOP : "<< ((double) 2 * n * n * k * i_max)/ss/1e9 <<endl;

  device.Memcpy(A_h,C,sa,MemcpyDeviceToHost);

  device.Free(A);
  device.Free(B);
  device.Free(C);

  free(A_h);
  free(B_h);
}
