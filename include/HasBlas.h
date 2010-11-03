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
