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

    void ssteqr_(char *compz, int &n, float *d, float *e, 
        float *z, int &ldz, float *work, int &info);

    void dsteqr_(char *compz, int &n, double *d, double *e, 
        double *z, int &ldz, double *work, int &info);
    
#ifdef __cplusplus
}
#endif

void steqr(char *compz, int &n, float *d, float *e, float *z, 
    int &ldz, float *work, int &info)
{
    ssteqr(compz, n, d, e, z, ldz, work, info);
}

void steqr(char *compz, int &n, double *d, double *e, double *z,
    int &ldz, double *work, int &info)
{
    dsteqr(compz, n, d, e, z, ldz, work, info);
}
