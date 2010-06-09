#ifndef _VESBLAS_H
#define _VESBLAS_H
   
#ifndef HAS_MKL_LIB
    
//#include <cblas.h>
#ifdef (__cplusplus)
extern "C"{
#endif

#define sgemm sgemm_
    extern "C"{
        void sgemm_(const char* TRANSA, const char* TRANSB,
            const int* M, const int* N, const int* K,
            const float* ALPHA,
            const float* A, const int* LDA, const float* B, const int* LDB, 
            const float* BETA, float* C, const int* LDC); 
    }
    
#ifdef (__cplusplus)
}
#endif

#else 

#include <mkl.h>
void steqr(char *compz, int &n, float *d, float *e, float *z, int &ldz, float *work, int &info)
{
    ssteqr(compz, &n, d, e, z, &ldz, work, &info);
}

void steqr(char *compz, int &n, double *d, double *e, double *z, int &ldz, double *work, int &info)
{
    dsteqr(compz, &n, d, e, z, &ldz, work, &info);
}

#endif //HAS_MKL_LIB
#endif //_VESBLAS_H_

