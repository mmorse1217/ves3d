#include "Device.h"

typedef float real;

int main(int argc, char **argv)
{
    COUT("\n ==============================\n"
        <<"  Blas Test:"
        <<"\n ==============================\n");
    sleep(1);

    int k = 2*2048;
    int n = 2*2048;

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
        Device<CPU> cpu(0);

        real *A = (real*) cpu.Malloc(sa * sizeof(real));
        real *B = (real*) cpu.Malloc(sb * sizeof(real));
        real *C = (real*) cpu.Malloc(sc * sizeof(real));

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

        real *A = (real*) gpu.Malloc(sa * sizeof(real));
        real *B = (real*) gpu.Malloc(sb * sizeof(real));
        real *C = (real*) gpu.Malloc(sc * sizeof(real));

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

    PROFILEREPORT(SortFlopRate);
    sleep(.5);
}
