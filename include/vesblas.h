#ifndef _VESBLAS_H
#define _VESBLAS_H



extern "C"{
#include <cblas.h>
}
//#define sgemm(TA,TB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc)							\
//	cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,*(M),*(N),*(K),	\
//							*(alpha),A,*(lda),B,*(ldb),*(beta),C,*(ldc))

#define sgemm sgemm_
extern "C"{
 void sgemm_(const char* TRANSA, const char* TRANSB,
						 const int* M, const int* N, const int* K,
						 const float* ALPHA,
						 const float* A, const int* LDA, const float* B, const int* LDB, 
             const float* BETA, float* C, const int* LDC); 
}

#endif
