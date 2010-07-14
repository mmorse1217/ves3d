#ifndef _VESBLAS_H
#define _VESBLAS_H
   
#ifndef HAS_MKL_LIB
#include <cblas.h>
// #include "f2c.h"
// #include "clapack.h"

#ifdef __cplusplus
extern "C"{
#endif

#define sgemm sgemm_
    extern "C"{
        void sgemm_(const char* TRANSA, const char* TRANSB,
            const int* M, const int* N, const int* K,
            const float* ALPHA, const float* A, const int* LDA, 
            const float* B, const int* LDB, 
            const float* BETA, float* C, const int* LDC); 
    }

#define dgemm dgemm_
    extern "C"{
        void dgemm_(const char* TRANSA, const char* TRANSB,
            const int* M, const int* N, const int* K,
            const double* ALPHA, const double* A, const int* LDA,
            const double* B, const int* LDB, 
            const double* BETA, double* C, const int* LDC); 
    }
    
// #define ssteqr ssteqr_
//     extern "C"{
//         void ssteqr_(char *compz, int &n, float *d, float *e, 
//             float *z, int &ldz, float *work, int &info);
//     }

// #define dsteqr dsteqr_
//     extern "C"{
//         void dsteqr_(char *compz, int &n, double *d, double *e, 
//             double *z, int &ldz, double *work, int &info);
//     }
    
#ifdef __cplusplus
}
#endif

#else 

#include <mkl.h>
#include <mkl_spblas.h>

void coomm(char *transa, int *m, int *n, int *k, float *alpha, 
    char *matdescra, float *val, int *rowind, int *colind, 
    int *nnz, float *b, int *ldb, float *beta, float *c, int *ldc)
{
    mkl_scoomm(transa, m, n, k, alpha, matdescra, val, rowind, 
        colind, nnz, b, ldb, beta, c, ldc);
}

void coomm(char *transa, int *m, int *n, int *k, double *alpha,
    char *matdescra, double *val, int *rowind, int *colind,
    int *nnz, double *b, int *ldb, double *beta, double *c, int *ldc)
{
    mkl_dcoomm(transa, m, n, k, alpha, matdescra, val, rowind, 
        colind, nnz, b, ldb, beta, c, ldc);
}

#endif //HAS_MKL_LIB

void steqr(char *compz, int &n, float *d, float *e, float *z, 
    int &ldz, float *work, int &info)
{
    ssteqr(compz, &n, d, e, z, &ldz, work, &info);
}

void steqr(char *compz, int &n, double *d, double *e, double *z,
    int &ldz, double *work, int &info)
{
    dsteqr(compz, &n, d, e, z, &ldz, work, &info);
}


#ifdef GPU_ACTIVE
#include "cublas.h"
#endif //GPU_ACTIVE

#endif //_VESBLAS_H_

