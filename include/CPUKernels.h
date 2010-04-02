#ifndef _CPUKERNELS_H_
#define _CPUKERNELS_H_

#include <emmintrin.h>
#include <omp.h>
#include <math.h>

void DirectStokesKernel(int stride, int n_surfs, int trg_idx_head,int trg_idx_tail, const float *qw, 
    const float *trg, const float *src, const float *den, float *pot);

void DirectStokesKernel_Noqw(int stride, int n_surfs, int trg_idx_head,int trg_idx_tail, 
    const float *trg, const float *src, const float *den, float *pot);

void DirectStokesSSE(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, 
    const float *qw, const float *trg, const float *src, const float *den, float *pot);

#endif//_CPUKERNELS_H_
