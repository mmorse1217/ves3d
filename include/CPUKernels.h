#ifndef _CPUKERNELS_H_
#define _CPUKERNELS_H_

#include <emmintrin.h>
#include <omp.h>
#include <math.h>
#include "Logger.h"
#include "PVFMMInterface.h"

///Single precision
void DirectStokesKernel(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const float *qw, const float *trg,
    const float *src, const float *den, float *pot);

void DirectStokesKernel_Noqw(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const float *trg, const float *src,
    const float *den, float *pot);

///Double precision
void DirectStokesKernel(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const double *qw, const double *trg,
    const double *src, const double *den, double *pot);

void DirectStokesKernel_Noqw(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const double *trg, const double *src,
    const double *den, double *pot);

///Single precision -- SSE instructions
void DirectStokesSSE(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const float *qw, const float *trg,
    const float *src, const float *den, float *pot);

void DirectStokesSSE_Noqw(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const float *trg, const float *src,
    const float *den, float *pot);

///Double precision -- SSE instructions
void DirectStokesSSE(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const double *qw, const double *trg,
    const double *src, const double *den, double *pot);

void DirectStokesSSE_Noqw(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const double *trg, const double *src,
    const double *den, double *pot);

///All to all interaction
void StokesAlltoAll(const float *src, const float *den, size_t np,
    float *pot, void**);
void StokesAlltoAll(const double *src, const double *den, size_t np,
    double *pot, void**);





void DirectStokesDoubleLayerKernel(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const float *qw, const float *trg,
    const float *src, const float *norm, const float *den, float *pot);

void DirectStokesDoubleLayerKernel(int src_stride, int trg_stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const double *qw, const double *trg,
    const double *src, const double *norm, const double *den, double *pot);

#endif//_CPUKERNELS_H_
