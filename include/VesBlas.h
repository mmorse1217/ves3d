#ifndef _VESBLAS_H
#define _VESBLAS_H
   
#ifdef HAS_MKL
  #include "HasMkl.h"

#elif  HAS_ATLAS
  #include "HasAtlas.h"

#elif  HAS_BLAS
  #include "HasBlas.h"

#endif

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

#endif //_VESBLAS_H_

