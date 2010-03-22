/**
 * @file   Device.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Feb 23 14:59:15 2010
 * 
 * @brief  The abstract class for the device class.
 */

#ifndef _DEVICE_H_
#define _DEVICE_H_

#include <sys/time.h>
#define get_seconds()   (gettimeofday(&tp, &tzp), \
    (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0)

struct timeval  tp;
struct timezone tzp;

//Forward declaration for the friend function
template<typename T> class Device;

// The comparison operator for the device class
template<typename Tlhs, typename Trhs>
inline bool operator==(const Device<Tlhs> &rhs, const Device<Trhs> &lhs)
{
    return( (void*) &rhs == (void*) &lhs); 
}

///The enum types for the memory copying action.
enum MemcpyKind {MemcpyHostToHost, MemcpyHostToDevice, MemcpyDeviceToHost, MemcpyDeviceToDevice};

///The enum type for the reordering of the points
enum CoordinateOrder {PointMajor, AxisMajor};
/**
 *  This class provides memory operations and also device specific
 *  optimized kernels <b>[Note the notational conventions in the
 *  Detailed Description section]</b>. The notational convention is
 *  that <tt>a</tt> and <tt>b</tt> are single scalars of type T;
 *  <tt>x</tt>, <tt>y</tt>, and <tt>z</tt> are scalar fields,
 *  i.e. arrays of scalar of type T with length
 *  <tt>stride*num_surfs</tt>; and <tt>u</tt>, <tt>v</tt>, and
 *  <tt>w</tt> are vector fields, i.e. arrays of scalar of type T with
 *  length <tt>3*stride*num_surfs</tt>. This convention is used
 *  throughout the class for both the naming of methods and naming the
 *  methods' parameters. All arrays are assumed to be residing on the
 *  device and there is no need of copying to the device etc.
 *
 *  @todo The tester does not cover the comparison operator, the
 *  ID'ing mechanism, and the sqrt
 */
template<typename T> class Device
{
  public:
    ///Memory allocation, since the type is defined as a template
    ///parameter, the size is different from the original malloc() as
    ///it <b>does not need to be multiplied by the sizeof(T)</b>.
    virtual T* Malloc(unsigned long int length) = 0;

    ///Freeing memory. 
    virtual void Free(T* ptr) = 0;

    ///Memory allocation and zero initialization for an array in
    ///memory, since the type is defined as a template parameter,
    ///there is no size argument.
    virtual T* Calloc(unsigned long int num) = 0;

    ///Copies the memory location form source to the destination. The
    ///kind of copy is from the enum type MemcpyKind and is <code>
    ///MemcpyHostToHost, MemcpyHostToDevice, MemcpyDeviceToHost, or
    ///MemcpyDeviceToDevice</code>.
    virtual T* Memcpy(T* destination, const T* source, unsigned long int num, enum MemcpyKind kind) = 0;
 
    ///Geometric dot product of two (Cartesian) vectors. 
    virtual T* DotProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* x_out) = 0;

    ///Geometric cross product of two (Cartesian) vectors.
    virtual T* CrossProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* w_out) = 0;

    ///Square root operator
    virtual T* Sqrt(const T* x_in, int stride, int num_surfs, T* sqrt_out) = 0;

    ///Element-wise inverse (of a scalar field).
    virtual T* xInv(const T* x_in, int stride, int num_surfs, T* xInv_out) = 0;

    ///Element-wise multiplication of scalar fields
    virtual T* xy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xy_out) = 0;

    ///Element-wise division of scalar fields.
    virtual T* xyInv(const T* x_in, const T* y_in, int stride, int num_surfs, T* xyInv_out) = 0;

    ///Element-wise scaling of a vector field by a scalar fields.
    virtual T* uyInv(const T* u_in, const T* y_in, int stride, int num_surfs, T* uyInv_out) = 0;
    
    ///Scaling and addition of an scalar field to another field.
    virtual T* axpy(T a_in, const T* x_in, const T* y_in, int stride, int num_surfs , T* axpy_out) = 0;
    
    ///Scaling of an array and addition of a single scalar. 
    virtual T* axpb(T a_in, const T*  x_in, T b_in, int stride, int num_surfs , T*  axpb_out) = 0;
    
    ///Element-wise scaling and addition.
    virtual T* avpw(const T* a_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  avpw_out) = 0;
 
    ///Element-wise scaling and addition.
    virtual T* xvpw(const T* x_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  xvpw_out) = 0;
    
    ///Element-wise scaling and addition.
    virtual T* xvpb(const T* x_in, const T*  v_in, T b_in, int stride, int num_surfs, T*  xvpb_out) = 0;

    ///Smooth integral--reduction
    virtual T* Reduce(const T *x_in, const T *w_in, const T *quad_w_in, int stride, int num_surfs, T  *int_x_dw) = 0;

    virtual T* gemm(const char *transA, const char *transB, const int *m, const int *n, const int *k, const T *alpha, 
		    const T *A, const int *lda, const T *B, const int *ldb, const T *beta, T *C, const int *ldc)=0;

    virtual T* CircShift(const T *arr_in, int n_vecs, int vec_length, int shift, T *arr_out) = 0;
    
    virtual void DirectStokes(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, 
        const T *qw, const T *trg, const T *src, const T *den, T *pot) = 0;

    virtual T* ShufflePoints(T *x_in, CoordinateOrder order_in, int stride, int n_surfs, T *x_out) = 0;

    ///SHT size
    virtual void InitializeSHT(int p, int p_up) = 0;
    
    ///SHT Analysis
    virtual void ShAna(const T *x_in, T *work_arr, int p, int num_funs, T *sht_out) = 0;
    
    ///SHT Synthesis
    virtual void ShSyn(const T *shc_in, T *work_arr, int p, int num_funs, T *y_out) = 0;
    
    ///SHT First derivative u
    virtual void ShSynDu(const T *shc_in, T *work_arr, int p, int num_funs, T *xu_out) = 0;
    
    ///SHT First derivative v
    virtual void ShSynDv(const T *shc_in, T *work_arr, int p, int num_funs, T *xv_out) = 0;
    
    ///SHT Second derivative uu
    virtual void ShSynDuu(const T *shc_in, T *work_arr, int p, int num_funs, T *xuu_out) = 0;
    
    ///SHT Second derivative vv
    virtual void ShSynDvv(const T *shc_in, T *work_arr, int p, int num_funs, T *xvv_out) = 0;
    
    ///SHT Second derivative uv
    virtual void ShSynDuv(const T *shc_in, T *work_arr, int p, int num_funs, T *xuv_out) = 0;

    virtual void AllDerivatives(const T *x_in, T *work_arr, int p, int num_funs, T* shc_x, T *Dux_out, T *Dvx_out, 
        T *Duux_out, T *Duvx_out, T *Dvvx_out) = 0;

    virtual void FirstDerivatives(const T *x_in, T *work_arr, int p, int num_funs, T* shc_x, T *Dux_out, T *Dvx_out) = 0;

    ///Filter
    virtual void Filter(int p, int n_funs, const T *x_in, const T *alpha, T* work_arr, T *shc_out, T *x_out) = 0;
    virtual void ScaleFreqs(int p, int n_funs, const T *inputs, const T *alphas, T *outputs) = 0;
    virtual void Resample(int p, int n_funs, int q, const T *shc_p, T *shc_q) = 0;

    ///The comparison operator for the device class
    template<typename Tlhs,typename Trhs>
    friend bool operator==(const Device<Tlhs> &rhs, const Device<Trhs> &lhs);
};

#endif //_DEVICE_H_
