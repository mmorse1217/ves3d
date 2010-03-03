/**
 * @file   DeviceCPU.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Feb 28 20:14:40 2010
 */

#ifndef _DEVICECPU_H_
#define _DEVICECPU_H_

#include "Device.h"
#include <iostream>
#include <cassert>
#include <math.h>

using namespace std;

///The CPU subclass of the Device class.
template<typename T>
class DeviceCPU : public Device<T>
{
  public:
    //memory operators
    virtual T* Malloc(unsigned long int length);
    virtual void  Free(T* ptr);
    virtual T* Calloc(unsigned long int num);

    ///@todo Memcpy is incomplete
    virtual T* Memcpy (T* destination, const T* source, unsigned long int num, enum MemcpyKind kind = MemcpyHostToHost);

    //Algebraic operators
    virtual T* DotProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* x_out);
    virtual T* CrossProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* w_out);

    virtual T* Sqrt(const T* x_in, int stride, int num_surfs, T* sqrt_out);
    virtual T* xInv(const T* x_in, int stride, int num_surfs, T* xInv_out);
    virtual T* xy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xy_out);
    virtual T* xyInv(const T* x_in, const T* y_in, int stride, int num_surfs, T* xyInv_out);

    virtual T* axpy(T a_in, const T* x_in, const T* y_in, int stride, int num_surfs , T* axpy_out);
    virtual T* axpb(T a_in, const T*  x_in, T b_in, int stride, int num_surfs , T*  axpb_out);
    virtual T* xvpw(const T* x_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  xvpw_out);
    virtual T* xvpb(const T* x_in, const T*  v_in, T b_in, int stride, int num_surfs, T*  xvpb_out);
};

#include "DeviceCPU.cc"
#endif //_DEVICECPU_H_
