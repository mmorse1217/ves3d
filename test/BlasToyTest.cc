#include <unistd.h>  //for sleep()
#include "Device.h"

typedef float real;

int main(int argc, char **argv)
{
    PROFILESTART();
    COUT("\n ==============================\n"
        <<"  Blas Test:"
        <<"\n ==============================\n");
    sleep(1);

    int k = 4*2048;
    int n = 4*2048;

    size_t sa = n*n;
    size_t sb = n*k;
    size_t sc = n*k;

    real *A_h = (real*) malloc(sa * sizeof(real));
    real *B_h = (real*) malloc(sb * sizeof(real));

    real alpha = 1.0;
    real beta = 0.0;

    for(int ii=0;ii<n;++ii)
        for(int jj=0;jj<n;++jj)
            A_h[ii*n+jj] = drand48();

    for(int ii=0;ii<n;++ii)
        for(int jj=0;jj<k;++jj)
            B_h[ii*k+jj] = drand48();

    {// cpu
      COUT(" - Testing CPU"<<endl);
        Device<CPU> cpu(0);

        real *A = (real*) cpu.Malloc(sa * sizeof(real));
        real *B = (real*) cpu.Malloc(sb * sizeof(real));
        real *C = (real*) cpu.Malloc(sc * sizeof(real));

        cpu.Memcpy(A,A_h,sa,Device<CPU>::MemcpyHostToDevice);
        cpu.Memcpy(B,B_h,sa,Device<CPU>::MemcpyHostToDevice);
        cpu.gemm("N", "N", &n, &k, &n, &alpha, A, &n, B, &n, &beta, C, &n);

        cpu.Free(A);
        cpu.Free(B);
        cpu.Free(C);
    }

#ifdef GPU_ACTIVE
    {
      COUT(" - Testing GPU"<<endl);
        Device<GPU> gpu(0);

        real *A = (real*) gpu.Malloc(sa * sizeof(real));
        real *B = (real*) gpu.Malloc(sb * sizeof(real));
        real *C = (real*) gpu.Malloc(sc * sizeof(real));

        gpu.Memcpy(A,A_h,sa,Device<GPU>::MemcpyHostToDevice);
        gpu.Memcpy(B,B_h,sa,Device<GPU>::MemcpyHostToDevice);
        gpu.gemm("N", "N", &n, &k, &n, &alpha, A, &n, B, &n, &beta, C, &n);

        gpu.Free(A);
        gpu.Free(B);
        gpu.Free(C);
    }
#endif //GPU_ACTIVE

    free(A_h);
    free(B_h);

    PROFILEEND("",0);
    PROFILEREPORT(SortFlop);
    sleep(.5);
}
