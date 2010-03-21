#ifndef _VESBLAS_H
#define _VESBLAS_H

#if defined (__cplusplus)
extern "C"{
#endif


#ifndef HAS_MKL_LIB

#include <cblas.h>
#define sgemm sgemm_
extern "C"{
 void sgemm_(const char* TRANSA, const char* TRANSB,
						 const int* M, const int* N, const int* K,
						 const float* ALPHA,
						 const float* A, const int* LDA, const float* B, const int* LDB, 
             const float* BETA, float* C, const int* LDC); 
}
#else 
#include <mkl.h>
#endif //HAS_MKL_LIB

#if defined (__cplusplus)
}
#endif
#endif
