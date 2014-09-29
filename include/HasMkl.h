#include <mkl.h>
#include <mkl_spblas.h>

#define NO_SPARSE_MATVEC false

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

