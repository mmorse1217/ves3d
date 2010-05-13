/**
 * @file   DeviceCPU.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Feb 28 20:14:40 2010
 */

#ifndef _DEVICECPU_H_
#define _DEVICECPU_H_

#include "DataIO.h"
#include "Device.h"
#include <iostream>
#include <cassert>
#include <cstring>
//#include <math.h>
//#include <memory.h>
#include "BlasSht.h"
#include "CPUKernels.h"
//#include <emmintrin.h>
#include <omp.h>
#include "OperatorsMats.h"

using namespace std;

struct CpuTime{
    static double gemm_time;
    static double stokes_time;
    static double xvpb_time;
    static double xy_time;
    static double DotProduct_time;
    static double Shift_time;
};

double CpuTime::gemm_time = 0;
double CpuTime::stokes_time = 0;
double CpuTime::xvpb_time = 0;
double CpuTime::xy_time = 0;
double CpuTime::DotProduct_time = 0;
double CpuTime::Shift_time = 0;
 
///The CPU subclass of the Device class.
template<typename T>
class DeviceCPU : public Device<T>
{
  public:
    //  BlasSht sht_;
    //BlasSht sht_up_sample_;
    //int p_, p_up_;

  public:
    DeviceCPU();
    ~DeviceCPU();

    //memory operators
    virtual T* Malloc(unsigned long int length) const;
    virtual void Free(T* ptr) const;
    virtual T* Calloc(unsigned long int num) const;

    ///@todo Memcpy is incomplete
    virtual T* Memcpy (T* destination, const T* source, unsigned long int num, enum MemcpyKind kind) const;
    virtual T* Memset (T *ptr, int value, size_t num) const;

    //Algebraic operators
    virtual T* DotProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* x_out) const;
    virtual T* CrossProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* w_out) const;

    virtual T* Sqrt(const T* x_in, int stride, int num_surfs, T* sqrt_out) const;
    virtual T* xInv(const T* x_in, int stride, int num_surfs, T* xInv_out) const;
    virtual T* xy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xy_out) const;
    virtual T* xyInv(const T* x_in, const T* y_in, int stride, int num_surfs, T* xyInv_out) const;
    virtual T* uyInv(const T* u_in, const T* y_in, int stride, int num_surfs, T* uyInv_out) const;

    virtual T* axpy(T a_in, const T* x_in, const T* y_in, int stride, int num_surfs , T* axpy_out) const;
    virtual T* axpb(T a_in, const T*  x_in, T b_in, int stride, int num_surfs , T*  axpb_out) const;
    virtual T* avpw(const T* a_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  avpw_out) const;
    virtual T* xvpw(const T* x_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  xvpw_out) const;
    virtual T* xvpb(const T* x_in, const T*  v_in, T b_in, int stride, int num_surfs, T*  xvpb_out) const;
    
    virtual T* Reduce(const T *x_in, const T *w_in, const T *quad_w_in, int stride, int num_surfs, T  *int_x_dw) const;

    virtual T* gemm(const char *transA, const char *transB, const int *m, const int *n, const int *k, const T *alpha, 
		    const T *A, const int *lda, const T *B, const int *ldb, const T *beta, T *C, const int *ldc) const;

    virtual T* CircShift(const T *arr_in, int n_vecs, int vec_length, int shift, T *arr_out) const;
    virtual void DirectStokes(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail,
        const T *qw, const T *trg, const T *src, const T *den, T *pot) const;
    
    virtual T* ShufflePoints(T *x_in, CoordinateOrder order_in, int stride, int n_surfs, T *x_out) const;

    virtual T Max(T *x_in, int length) const;
    //SHT
//     virtual void InitializeSHT(OperatorsMats<T> &mats);

//     virtual void ShAna(   const T *x_in  , T *work_arr, int p, int num_funs, T *shc_out) const;
//     virtual void ShSyn(   const T *shc_in, T *work_arr, int p, int num_funs, T *x_out  ) const;
//     virtual void ShSynDu( const T *shc_in, T *work_arr, int p, int num_funs, T *xu_out ) const;
//     virtual void ShSynDv( const T *shc_in, T *work_arr, int p, int num_funs, T *xv_out ) const;
//     virtual void ShSynDuu(const T *shc_in, T *work_arr, int p, int num_funs, T *xuu_out) const;
//     virtual void ShSynDvv(const T *shc_in, T *work_arr, int p, int num_funs, T *xvv_out) const;
//     virtual void ShSynDuv(const T *shc_in, T *work_arr, int p, int num_funs, T *xuv_out) const;

//     virtual void AllDerivatives  (const T *x_in, T *work_arr, int p, int num_funs, 
//         T *shc_x, T *Dux_out, T *Dvx_out,T *Duux_out, T *Duvx_out, T *Dvvx_out) const;
//     virtual void FirstDerivatives(const T *x_in, T *work_arr, int p, int num_funs, T *shc_x, T *Dux_out, T *Dvx_out) const;

//     virtual void Filter(int p, int n_funs, const T *x_in, const T *alpha, T* work_arr, T *shc_out, T *x_out) const;
//     virtual void ScaleFreqs(int p, int n_funs, const T *inputs, const T *alphas, T *outputs) const;
//     virtual void Resample(int p, int n_funs, int q, const T *shc_p, T *shc_q) const;
//     virtual void InterpSh(int p, int n_funs, const T *x_in, T* work_arr, T *shc, int q, T *x_out) const;
};

#include "DeviceCPU.cc"
#endif //_DEVICECPU_H_
