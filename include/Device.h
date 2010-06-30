/**
 * @file   Device.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Feb 23 14:59:15 2010
 */

#ifndef _DEVICE_H_
#define _DEVICE_H_

#include <cstring> //for the memcpy and memset
#include <cassert>
#include <iostream> 
#include <cmath>
#include <omp.h>
#include "VesBlas.h"
#include "Logger.h"
#include "CPUKernels.h"
#include "enums.h"

#ifdef GPU_ACTIVE
#include "cuda_runtime.h"
#include "CudaKernels.h"
#endif //GPU_ACTIVE

///
#define EMPTY_GRID make_pair(-1,-1)
#define DIM 3

// Forward declaration
template<enum DeviceType DT> class Device;

/// The comparison operator for the device class 
template<enum DeviceType DTlhs, enum DeviceType DTrhs>                         
inline bool operator==(const Device<DTlhs> &lhs, const Device<DTrhs> &rhs);

///
inline pair<int, int> gridDimOf(int sh_order);

/**
 * This class provides memory operations and device specific optimized
 * kernels <b>[Note the notational conventions in the Detailed
 * Description section]</b>. The notational convention is that
 * <tt>a</tt> and <tt>b</tt> are single scalars of type T; <tt>x</tt>,
 * <tt>y</tt>, and <tt>z</tt> are scalar fields, i.e. arrays of scalar
 * of type T with length <tt>stride*num_vecs</tt>; and <tt>u</tt>,
 * <tt>v</tt>, and <tt>w</tt> are vector fields, i.e. arrays of scalar
 * of type T with length <tt>3*stride*num_vecs</tt>. This convention
 * is used throughout the class for both the naming of methods and
 * naming the methods' parameters. All arrays are assumed to be
 * residing on the device.
 *
 * The behavior of <tt>Device<GPU></tt> depends on the variable
 * <tt>GPU_ACTIVE</tt>. If <tt>GPU_ACTIVE</tt> is defined, the class
 * calls the GPU specialized methods, otherwise the generic methods
 * (CPU) are called.
 */
template<enum DeviceType DT = CPU>
class Device
{
  public:
    ///The constructor of the class. device_id is a user-specified ID
    ///and err is the return error type corresponding to the
    ///instantiation action. Both of these parameter are optional.
    explicit Device(int device_id = 0, enum DeviceError *err = NULL);

    ///Destructor
    ~Device();

    ///Memory allocation.
    void* Malloc(size_t length) const;

    ///Freeing memory. 
    void Free(void* ptr) const;

    ///Memory allocation and zero initialization for an array in
    ///memory.
    void* Calloc(size_t num, size_t size) const;

    ///Copies the memory location form source to the destination. The
    ///kind of copy is from the enum type MemcpyKind that includes
    ///<code> MemcpyHostToHost, MemcpyHostToDevice,
    ///MemcpyDeviceToHost, or MemcpyDeviceToDevice</code>.
    void* Memcpy(void* destination, const void* source, 
        size_t num, enum MemcpyKind kind) const;
    
    ///Setting the memory location specified by ptr to the given value
    void* Memset(void *ptr, int value, size_t num) const;
    
    ///Geometric dot product of two (Cartesian) vectors. 
    template<typename T>
    T* DotProduct(const T* u_in, const T* v_in, size_t stride, 
        size_t n_vecs, T* x_out) const;

    ///Geometric cross product of two (Cartesian) vectors.
    template<typename T>
    T* CrossProduct(const T* u_in, const T* v_in, size_t stride, 
        size_t n_vecs, T* w_out) const;

    ///Square root operator.
    template<typename T>
    T* Sqrt(const T* x_in, size_t length, T* sqrt_out) const;

    template<typename T>
    T* ax(const T* a, const T* x, size_t length, size_t n_vecs, T* ax_out) const;

    ///Element-wise multiplication of scalar fields.
    template<typename T>
    T* xy(const T* x_in, const T* y_in, size_t length, T* xy_out) const;

    ///Element-wise division of scalar fields.
    template<typename T>
    T* xyInv(const T* x_in, const T* y_in, size_t length, 
        T* xyInv_out) const;

    ///Element-wise scaling of a vector field by a scalar fields.
    template<typename T>
    T* uyInv(const T* u_in, const T* y_in, size_t stride, size_t n_vecs, 
        T* uyInv_out) const;
    
    ///Scaling and addition of an scalar field to another field.
    template<typename T>
    T* axpy(T a_in, const T* x_in, const T* y_in, size_t length, 
        T* axpy_out) const;
    
    template<typename T>
    T* apx(T* a_in, const T* x_in, size_t stride, 
        size_t n_subs, T* axpy_out) const;
    
    ///Element-wise scaling and addition.
    template<typename T>
    T* avpw(const T* a_in, const T* v_in, const T* w_in, 
        size_t stride, size_t n_vecs, T*  avpw_out) const;
 
    ///Element-wise scaling and addition.
    template<typename T>
    T* xvpw(const T* x_in, const T*  v_in, const T*  w_in, 
        size_t stride, size_t n_vecs, T*  xvpw_out) const;

    ///Smooth integral (reduction) for multidimensional fields.
    template<typename T>
    T* Reduce(const T *x_in, const int x_dim, const T *w_in, const T *quad_w_in, 
        const size_t stride, const size_t ns, T *x_dw) const;
    
    ///General matrix-matrix multiplication. consult BLAS
    ///documentation for the detail of the syntax.
    float* gemm(const char *transA, const char *transB, const int *m, 
        const int *n, const int *k, const float *alpha, const float *A, 
        const int *lda, const float *B, const int *ldb, const float *beta, 
        float *C, const int *ldc) const;

    ///General matrix-matrix multiplication. consult BLAS
    ///documentation for the detail of the syntax.
    double* gemm(const char *transA, const char *transB, const int *m, 
        const int *n, const int *k, const double *alpha, const double *A, 
        const int *lda, const double *B, const int *ldb, const double *beta, 
        double *C, const int *ldc) const;
    
    ///Direct stokes integration.
    void DirectStokes(const float *src, const float *den, const float *qw, 
        size_t stride, size_t n_surfs, const float *trg, size_t trg_idx_head, 
        size_t trg_idx_tail, float *pot) const;

    ///Direct stokes integration.
    void DirectStokes(const double *src, const double *den, const double *qw, 
        size_t stride, size_t n_surfs, const double *trg, size_t trg_idx_head, 
        size_t trg_idx_tail, double *pot) const;

    ///The max of an array.
    template<typename T>
    T MaxAbs(T *x_in, size_t length) const;

    ///Transpose of a matrix
    template<typename T>
    T* Transpose(const T *in, size_t height, size_t width, T *out) const;

    template<typename T>
    T AlgebraicDot(const T* x, const T* y, size_t length) const;

    //The comparison operator 
    template<enum DeviceType DTlhs, enum DeviceType DTrhs>                         
    friend bool operator==(const Device<DTlhs> &lhs, const Device<DTrhs> &rhs);

  private:
    ///The declaration of copy constructor. It is private to avoid
    ///pass by value. There is no implementation for this method.
    Device(const Device<DT> &device_in);
    
    /* The declaration of the assignment operator, it is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    Device<DT>& operator=(const Device<DT> &device_in);
};

#include "DeviceCPU.cc"

#ifdef GPU_ACTIVE
#include "DeviceGPU.cc"
#endif

#endif //_DEVICE_H_
