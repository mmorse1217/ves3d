#ifndef _BLAS_H_
#define _BLAS_H_

/*
#define sgemm sgemm_
#define dgemm dgemm_   
#define ssteqr ssteqr_
#define dsteqr dsteqr_

#ifdef __cplusplus
extern "C"{
#endif
    //#include <cblas.h>
    //#include <clapack.h>

    void sgemm_(const char* TRANSA, const char* TRANSB,
        const int* M, const int* N, const int* K,
        const float* ALPHA, const float* A, const int* LDA, 
        const float* B, const int* LDB, 
        const float* BETA, float* C, const int* LDC);  

    void dgemm_(const char* TRANSA, const char* TRANSB,
        const int* M, const int* N, const int* K,
        const double* ALPHA, const double* A, const int* LDA,
        const double* B, const int* LDB, 
        const double* BETA, double* C, const int* LDC); 

    void ssteqr_(char *compz, const int *n, float *d, float *e, 
        float *z, const int *ldz, float *work, const int *info);

    void dsteqr_(char *compz, const int *n, double *d, double *e, 
        double *z, const int *ldz, double *work, const int *info);
    
#ifdef __cplusplus
}
#endif
*/

#include "mkl.h"

void Gemm(const char* transa, const char* transb,
    const int* m, const int* n, const int* k, const float* alpha, 
    const float* a, const int* lda, const float* b, const int* ldb, 
    const float* beta, float* c, const int* ldc)
{
  sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);  
}

void Gemm(const char* transa, const char* transb,
    const int* m, const int* n, const int* k, const double* alpha, 
    const double* a, const int* lda, const double* b, const int* ldb, 
    const double* beta, double* c, const int* ldc)
{
  dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);  
}

void Steqr(char *compz, int &n, float *d, float *e, float *z, 
    int &ldz, float *work, int &info)
{
  ssteqr(compz, &n, d, e, z, &ldz, work, &info);
}

void Steqr(char *compz, int &n, double *d, double *e, double *z,
    int &ldz, double *work, int &info)
{
  dsteqr(compz, &n, d, e, z, &ldz, work, &info);
}

#ifdef GPU_ACTIVE

#include "cublas.h"

void cugemm(const char *transa, const char *transb, 
    const int *m, const int *n, const int *k, const float *alpha, 
    const float *A, const int *lda, const float *B, const int *ldb, 
    const float *beta, float *C, const int *ldc)
{
    cublasSgemm(*transa, *transb, *m, *n, *k, *alpha, A, *lda, B, 
        *ldb, *beta, C, *ldc); 
}    

void cugemm(const char *transa, const char *transb, 
    const int *m, const int *n, const int *k, const double *alpha, 
    const double *A, const int *lda, const double *B, const int *ldb, 
    const double *beta, double *C, const int *ldc)
{
    cublasDgemm(*transa, *transb, *m, *n, *k, *alpha, A, *lda, B, 
        *ldb, *beta, C, *ldc); 
}    

#endif //GPU_ACTIVE

#endif //_BLAS_H_

