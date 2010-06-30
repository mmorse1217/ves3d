#include "Device.h"

int main(int argc, char **argv)
{
    COUT("\n ==============================\n"
        <<"  Blas Test:"
        <<"\n ==============================\n");
    sleep(1);
    
  int k = 4096;
  int n = 4096;

  size_t sa = n*n;
  size_t sb = n*k;
  size_t sc = n*k;

  float *A_h = (float*) malloc(sa * sizeof(float));
  float *B_h = (float*) malloc(sb * sizeof(float));

  float alpha = 1.0;
  float beta = 0.0;

  for(int ii=0;ii<n;++ii)
    for(int jj=0;jj<n;++jj)
      A_h[ii*n+jj] = drand48();

  for(int ii=0;ii<n;++ii)
    for(int jj=0;jj<k;++jj)
      B_h[ii*k+jj] = drand48();

  {// cpu
      Device<CPU> cpu(0);
      
      float *A = (float*) cpu.Malloc(sa * sizeof(float));
      float *B = (float*) cpu.Malloc(sb * sizeof(float));
      float *C = (float*) cpu.Malloc(sc * sizeof(float));
      
      cpu.Memcpy(A,A_h,sa,MemcpyHostToDevice);
      cpu.Memcpy(B,B_h,sa,MemcpyHostToDevice);
      cpu.gemm("N", "N", &n, &k, &n, &alpha, A, &n, B, &n, &beta, C, &n);
      
      cpu.Free(A);
      cpu.Free(B);
      cpu.Free(C);
  }

#ifdef GPU_ACTIVE
  {
      Device<GPU> gpu(0);

      float *A = (float*) gpu.Malloc(sa * sizeof(float));
      float *B = (float*) gpu.Malloc(sb * sizeof(float));
      float *C = (float*) gpu.Malloc(sc * sizeof(float));
      
      gpu.Memcpy(A,A_h,sa,MemcpyHostToDevice);
      gpu.Memcpy(B,B_h,sa,MemcpyHostToDevice);
      gpu.gemm("N", "N", &n, &k, &n, &alpha, A, &n, B, &n, &beta, C, &n);
      
      gpu.Free(A);
      gpu.Free(B);
      gpu.Free(C);
  }
#endif //GPU_ACTIVE

  free(A_h);
  free(B_h);

  PROFILEREPORT(SortTime);
}
