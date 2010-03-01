/**
 * @file   Device.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Feb 23 14:59:15 2010
 * 
 * @brief  The abstract class for the device class.
 */

#ifndef _DEVICE_H_
#define _DEVICE_H_

#include <iostream>
#include <cstring>
#include <cassert>
// cuda.h is included inside DeviceGPU implementation.

///The enum types for the memory copying action.
enum MemcpyKind {MemcpyHostToHost, MemcpyHostToDevice, MemcpyDeviceToHost, MemcpyDeviceToDevice};

template<typename T> class Device
{
  public:
    ///Memory allocation.
    virtual T* Malloc(size_t size) = 0;

//     ///Freeing memory. 
//     virtual void Free(T* ptr) = 0;

//     ///Memory allocation for an array in memory.
//     virtual T* Calloc(size_t num, size_t size) = 0;

//     ///Copies the memory location form source to the destination. The
//     ///kind of copy is from the enum type MemcpyKind and is
//     ///MemcpyHostToHost, MemcpyHostToDevice, MemcpyDeviceToHost, or
//     ///MemcpyDeviceToDevice.
//     virtual T* Memcpy (T* destination, const T* source, size_t num, enum MemcpyKind kind) = 0;

//     ///Geometric dot product of two (Cartesian) vectors residing on
//     ///the device, single precision.
//     virtual T* DotProduct(const T*  a_in, const T*  b_in, int stride, int num_surfs, T*  aDb_out) = 0;

//     ///Geometric dot product of two (Cartesian) vectors residing on
//     ///the device, double precision.
//     //virtual double* DotProduct(const double* a_in, const double* b_in, int stride, int num_surfs, double* aDb_out) = 0;

//     ///Geometric cross product of two (Cartesian) vectors residing on
//     ///the device, single precision.
//     virtual T*  CrossProduct(const T*  a_in, const T*  b_in, int stride, int num_surfs, T*  aCb_out) = 0;

//     ///Geometric cross product of two (Cartesian) vectors residing on
//     ///the device, double precision.
//     //virtual double* CrossProduct(const double* a_in, const double* b_in, int stride, int num_surfs, double* aCb_out) = 0;

//     ///Element-wise inverse (of a scalar field)
//     virtual T* xInv(const T* x_in, int stride, int num_surfs, T* xInv_out) = 0;
    
//     ///Scaling of an array and addition of scalar fields, single precision
//     virtual T* AxPy(T a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs , T*  axpy_out) = 0;

//     ///Scaling of an array and addition of scalar fields, single precision
//     virtual T* AxPy(T a_in, const T*  x_in, T y_in, int stride, int num_surfs , T*  axpy_out) = 0;
    
//     ///Element-wise scaling and addition, x and y vector field, a scalar field, single precision
//     virtual T* AxPy(const T* a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs, T*  axpy_out) = 0;
    
//     ///Element-wise scaling and addition, single precision
//     virtual T* AxPy(const T* a_in, const T*  x_in, T y_in, int stride, int num_surfs, T*  axpy_out) = 0;
   
//     ///Element-wise multiplication of scalar fields
//     virtual T* xTy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xTy_out) = 0;
};

///The GPU subclass of the Device class.
// class DeviceGPU : public Device
// {
//   public:
//     //memory operators
//     virtual void* Malloc(size_t size);
//     virtual void  Free(void* ptr);
//     virtual void* Calloc(size_t num, size_t size);
//     virtual void* Memcpy (void* destination, const void* source, size_t num, enum MemcpyKind kind);

//     //Algebraic operators
//     virtual T*  DotProduct(const T*  a_in, const T*  b_in, int stride, int num_surfs, T*  aDb_out);
//     virtual double* DotProduct(const double* a_in, const double* b_in, int stride, int num_surfs, double* aDb_out);

//     virtual T*  CrossProduct(const T*  a_in, const T*  b_in, int stride, int num_surfs, T*  aCb_out);
//     virtual double* CrossProduct(const double* a_in, const double* b_in, int stride, int num_surfs, double* aCb_out);

//     virtual T*  AxPy(T a_in, const T*  x_in, const T*  y_in, int stride, int num_scs , T*  axpy_out) = 0;
//     virtual T*  AxPy(T a_in, const T*  x_in, T y_in, int stride, int num_scs , T*  axpy_out) = 0;
//     virtual T*  AxPy(const T* a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs, T*  axpy_out) = 0;
//     virtual T*  AxPy(const T* a_in, const T*  x_in, T y_in, int stride, int num_surfs, T*  axpy_out) = 0;
   
    
//     virtual T* xTy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xTy_out) = 0;

//     virtual T* xDy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xDy_out) = 0;
// };

///The CPU subclass of the Device class.
// template<typename T>
// class DeviceCPU : public Device<T>
// {
//   public:
//     //memory operators
//     virtual T* Malloc(size_t size);
//     virtual void  Free(void* ptr);
//     virtual void* Calloc(size_t num, size_t size);
//     virtual void* Memcpy (void* destination, const void* source, size_t num, enum MemcpyKind kind = MemcpyHostToHost);

//     //Algebraic operators
//     virtual T* DotProduct(const T*  a_in, const T*  b_in, int stride, int num_surfs, T*  aDb_out);
//     //virtual T*  DotProduct(const T*  a_in, const T*  b_in, int stride, int num_surfs, T*  aDb_out);
//     //virtual double* DotProduct(const double* a_in, const double* b_in, int stride, int num_surfs, double* aDb_out);

//     virtual T*  CrossProduct(const T*  a_in, const T*  b_in, int stride, int num_surfs, T*  aCb_out);
//     //virtual double* CrossProduct(const double* a_in, const double* b_in, int stride, int num_surfs, double* aCb_out);

//     virtual T* xInv(const T* x_in, int stride, int num_surfs, T* xInv_out);

//     virtual T* AxPy(T a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs , T*  axpy_out);
//     virtual T* AxPy(T a_in, const T*  x_in, T y_in, int stride, int num_surfs , T*  axpy_out);
//     virtual T* AxPy(const T* a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs, T*  axpy_out);
//     virtual T* AxPy(const T* a_in, const T*  x_in, T y_in, int stride, int num_surfs, T*  axpy_out);
   
//     virtual T* xTy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xTy_out);
// };

// #include "DeviceCPU.cc"
#endif //_DEVICE_H_
