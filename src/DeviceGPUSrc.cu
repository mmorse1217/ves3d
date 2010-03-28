/**
 * @file   DeviceGPU.cu
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Mon Mar  1 13:09:06 2010
 */

template <typename T>
DeviceGPU<T>::DeviceGPU()
{
    ///@bug this should be unique
    cudaSetDevice(0);
}

template <typename T>
T* DeviceGPU<T>::Malloc(unsigned long int length)
{
    T* ptr = 0;
    cudaMalloc(&ptr, length * sizeof(T));
    return(ptr);
}

template <typename T>
void DeviceGPU<T>::Free(T* ptr)
{
    cudaFree(ptr);
}

template <typename T>
T* DeviceGPU<T>::Calloc(unsigned long int num)
{
    T* ptr = 0;
    cudaMalloc(&ptr, num * sizeof(T));
    cudaMemset(ptr, 0, num * sizeof(T));
    return ptr;
}

template <typename T>
T* DeviceGPU<T>::Memcpy (T* destination, const T* source, unsigned long int num, enum MemcpyKind kind)
{
    cudaMemcpyKind cuda_kind = static_cast<cudaMemcpyKind>(kind);
    cudaMemcpy(destination, source, sizeof(T) * num, cuda_kind);
    return destination;
}

#define DIM 3
#define BLOCK_HEIGHT 128

template <typename T> 
T* DeviceGPU<T>::DotProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* x_out)
{
    int GridDim = num_surfs;
    dotProdKernel<<<GridDim, BLOCK_HEIGHT>>>
        (u_in, v_in, stride, num_surfs, x_out);

    return x_out;
}

template <typename T> 
T* DeviceGPU<T>::CrossProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* w_out)
{
    return w_out;
}

template <typename T> 
T* DeviceGPU<T>::Sqrt(const T* x_in, int stride, int num_surfs, T* sqrt_out)
{
    return sqrt_out;
}

template <typename T> 
T* DeviceGPU<T>::xInv(const T* x_in, int stride, int num_surfs, T* xInv_out)
{
    return xInv_out;
}

template <typename T> 
T* DeviceGPU<T>::xy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xy_out)
{
    return xy_out;
}

template <typename T> 
T* DeviceGPU<T>::xyInv(const T* x_in, const T* y_in, int stride, int num_surfs, T* xyInv_out)
{
    return xyInv_out;
}

template <typename T>
T* DeviceGPU<T>::uyInv(const T* u_in, const T* y_in, int stride, int num_surfs, T* uyInv_out)
{
    return uyInv_out;
}

template <typename T> 
T* DeviceGPU<T>::axpy(T a_in, const T* x_in, const T* y_in, int stride, int num_surfs , T* axpy_out)
{
    return axpy_out;
}

template <typename T> 
T* DeviceGPU<T>::axpb(T a_in, const T*  x_in, T b_in, int stride, int num_surfs , T*  axpb_out)
{
    return axpb_out;
}

template <typename T> 
T* DeviceGPU<T>::avpw(const T* a_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  avpw_out)
{
    return avpw_out;
}

template <typename T> 
T* DeviceGPU<T>::xvpw(const T* x_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  xvpw_out)
{
    return xvpw_out;
}

template <typename T> 
T* DeviceGPU<T>::xvpb(const T* x_in, const T*  v_in, T b_in, int stride, int num_surfs, T*  xvpb_out)
{
    return xvpb_out;
}

template <typename T> 
T* DeviceGPU<T>::Reduce(const T *x_in, const T *w_in, const T *quad_w_in, int stride, int num_surfs, T  *int_x_dw)
{
    return int_x_dw;
}

template <typename T> 
T* DeviceGPU<T>::gemm(const char *transA, const char *transB, const int *m, const int *n, const int *k, const T *alpha, 
    const T *A, const int *lda, const T *B, const int *ldb, const T *beta, T *C, const int *ldc)
{
    return C;
}

template <typename T> 
T* DeviceGPU<T>::CircShift(const T *arr_in, int n_vecs, int vec_length, int shift, T *arr_out)
{
    return arr_out;
}


template <typename T>
void DeviceGPU<T>::DirectStokes(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, 
    const T *qw, const T *trg, const T *src, const T *den, T *pot)
{
}

template <typename T> 
T* DeviceGPU<T>::ShufflePoints(T *x_in, CoordinateOrder order_in, int stride, int n_surfs, T *x_out)
{
    return x_out;
}

template <typename T>
void DeviceGPU<T>::InitializeSHT(int p, int p_up)
{
}

template <typename T>
void DeviceGPU<T>::ShAna(const T *x_in, T *work_arr, int p, int num_funs, T *sht_out)
{
}

template <typename T>
void DeviceGPU<T>::ShSyn(const T *shc_in, T *work_arr, int p, int num_funs, T *y_out)
{
}

template <typename T>
void DeviceGPU<T>::ShSynDu(const T *shc_in, T *work_arr, int p, int num_funs, T *xu_out)
{
}

template <typename T>
void DeviceGPU<T>::ShSynDv(const T *shc_in, T *work_arr, int p, int num_funs, T *xv_out)
{
}

template <typename T>
void DeviceGPU<T>::ShSynDuu(const T *shc_in, T *work_arr, int p, int num_funs, T *xuu_out)
{
}

template <typename T>
void DeviceGPU<T>::ShSynDvv(const T *shc_in, T *work_arr, int p, int num_funs, T *xvv_out)
{
}

template <typename T>
void DeviceGPU<T>::ShSynDuv(const T *shc_in, T *work_arr, int p, int num_funs, T *xuv_out)
{
}

template <typename T>
void DeviceGPU<T>::AllDerivatives(const T *x_in, T *work_arr, int p, int num_funs, T* shc_x, T *Dux_out, T *Dvx_out, 
    T *Duux_out, T *Duvx_out, T *Dvvx_out)
{
}

template <typename T>
void DeviceGPU<T>::FirstDerivatives(const T *x_in, T *work_arr, int p, int num_funs, T* shc_x, T *Dux_out, T *Dvx_out)
{
}

template <typename T>
void DeviceGPU<T>::Filter(int p, int n_funs, const T *x_in, const T *alpha, T* work_arr, T *shc_out, T *x_out)
{
}

template <typename T>
void DeviceGPU<T>::ScaleFreqs(int p, int n_funs, const T *inputs, const T *alphas, T *outputs)
{
}

template <typename T>
void DeviceGPU<T>::Resample(int p, int n_funs, int q, const T *shc_p, T *shc_q)
{
}


