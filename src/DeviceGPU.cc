template<>
Device<GPU>::Device(int device_id, enum DeviceError *err)
{
    if(err != 0)
        *err = static_cast<DeviceError>(cudaSetDevice(device_id));
    else
        cudaSetDevice(device_id);

    cublasInit();
}

template<>
void* Device<GPU>::Malloc(size_t length) const
{
    PROFILESTART();
    void* ptr = 0;
    cudaMalloc(&ptr, length);
    PROFILEEND("GPU",0);
    return(ptr);
}
template<>
void Device<GPU>::Free(void* ptr) const
{
    PROFILESTART();
    cudaFree(ptr);
    PROFILEEND("GPU",0);
}
 
template<>
void* Device<GPU>::Calloc(size_t num, size_t size) const
{
    PROFILESTART();
    void* ptr = 0;
    cudaMalloc(&ptr, num * size);
    cudaMemset(ptr, 0, num * size);
    PROFILEEND("GPU",0);
    return ptr;
}

template<>
void* Device<GPU>::Memcpy(void* destination, const void* source, 
        size_t num, enum MemcpyKind kind) const
{
    PROFILESTART();
    cudaMemcpyKind cuda_kind = static_cast<cudaMemcpyKind>(kind);
    cudaMemcpy(destination, source, num, cuda_kind);
    PROFILEEND("GPU",0);
    return destination;
}
    
template<>
void* Device<GPU>::Memset(void *ptr, int value, size_t num) const
{
    PROFILESTART();
    cudaMemset(ptr, value, num);
    PROFILEEND("GPU",0);
    return(ptr);
}

template<>
template<typename T>  
T* Device<GPU>::DotProduct(const T* u_in, const T* v_in, size_t stride, 
    size_t n_vecs, T* x_out) const
{
    PROFILESTART();
    DotProductGpu(u_in, v_in, stride, n_vecs, x_out);
    PROFILEEND("GPU",0);
    return x_out;
}

template<>
template<typename T>  
T* Device<GPU>::CrossProduct(const T* u_in, const T* v_in, size_t stride, size_t num_surfs, T* w_out) const
{
    PROFILESTART();
    assert(DIM==3);
    CrossProductGpu(u_in, v_in, stride, num_surfs, w_out); 
    PROFILEEND("GPU",0);
    return w_out;
}

template<>
template<typename T>
T* Device<GPU>::Sqrt(const T* x_in, size_t length, T* sqrt_out) const
{
    PROFILESTART();
    SqrtGpu(x_in, length, sqrt_out);
    PROFILEEND("GPU",0);
    return sqrt_out;
}

template<>
template<typename T>
T* Device<GPU>::xy(const T* x_in, const T* y_in, size_t length, T* xy_out) const
{
    PROFILESTART();
    xyGpu(x_in, y_in, length, xy_out);
    PROFILEEND("GPU",0);
    return xy_out;
}

template<>
template<typename T>
T* Device<GPU>::xyInv(const T* x_in, const T* y_in, size_t length, T* xyInv_out) const
{
    PROFILESTART();
    if(x_in==NULL)
        InvGpu(y_in, length, xyInv_out);
    else
        xyInvGpu(x_in, y_in, length, xyInv_out);

    PROFILEEND("GPU",0);
    return xyInv_out;
}

template<>
template<typename T>
T*  Device<GPU>::uyInv(const T* u_in, const T* y_in, size_t stride, size_t num_surfs, T* uyInv_out) const
{
    PROFILESTART();
    assert(u!=NULL);
    uyInvGpu(u_in, y_in, stride, num_surfs, uyInv_out);
    PROFILEEND("GPU",0);
    return uyInv_out;
}

template<>
template<typename T>
T*  Device<GPU>::axpy(T a_in, const T*  x_in, const T*  y_in, size_t length, T*  axpy_out) const
{
    PROFILESTART();
    assert(x_in != NULL);
    if(y_in !=NULL)
        axpyGpu(a_in, x_in, y_in, length, axpy_out);
    else
        axpbGpu(a_in, x_in, (T) 0.0, length, axpy_out);
        
    PROFILEEND("GPU",0);
    return axpy_out;
}
template<>
template<typename T>
T*  Device<GPU>::avpw(const T* a_in, const T*  v_in, const T*  w_in, size_t stride, size_t num_surfs, T*  avpw_out) const
{
    PROFILESTART();
    if(w_in !=NULL)
        avpwGpu(a_in, v_in, w_in, stride, num_surfs, avpw_out);
    else
    {
        std::cerr<<"This kernel is not implemented for GPU"<<std::endl;
        abort();
    }

    PROFILEEND("GPU",0);
    return avpw_out;
}

template<>
template<typename T>
T*  Device<GPU>::xvpw(const T* x_in, const T*  v_in, const T*  w_in, size_t stride, size_t num_surfs, T*  xvpw_out) const
{
    PROFILESTART();
    if(w_in !=NULL)
        xvpwGpu(v_in, x_in, w_in, stride, num_surfs, xvpw_out);
    else
        xvpbGpu(v_in, x_in, (T) 0.0, stride, num_surfs, xvpw_out);
    
    PROFILEEND("GPU",0);
    return xvpw_out;
}

template<>
template<typename T>
T*  Device<GPU>::Reduce(const T *x_in, const T *w_in, const T *quad_w_in, size_t stride, size_t num_surfs, T  *int_x_dw) const
{
    PROFILESTART();
    ReduceGpu(x_in, w_in, quad_w_in, stride, num_surfs, int_x_dw);
    PROFILEEND("GPU",0);
    return int_x_dw;
}

template<>
float* Device<GPU>::gemm(const char *transA, const char *transB, 
    const int *m, const int *n, const int *k, const float *alpha, 
    const float *A, const int *lda, const float *B, const int *ldb, 
    const float *beta, float *C, const int *ldc) const
{
    PROFILESTART();
    cublasSgemm(*transA, *transB, *m, *n, *k, *alpha, A, *lda, B, *ldb, *beta, C, *ldc); 
    cudaThreadSynchronize();
    PROFILEEND("GPUs",(double) 2* (*k) * (*n) * (*m) + *(beta) * (*n) * (*m));
    return C;
}

template<>
double* Device<GPU>::gemm(const char *transA, const char *transB, 
    const int *m, const int *n, const int *k, const double *alpha, 
    const double *A, const int *lda, const double *B, const int *ldb, 
    const double *beta, double *C, const int *ldc) const
{
    std::cerr<<"The general matrix multiply is not implemented for the data type of double"<<std::endl;
    abort();
}

template<>
void Device<GPU>::DirectStokes(const float *src, const float *den, const float *qw, 
    size_t stride, size_t n_surfs, const float *trg, size_t trg_idx_head, 
    size_t trg_idx_tail, float *pot)
{
    PROFILESTART();
    cuda_stokes(stride, n_surfs, trg_idx_head, trg_idx_tail, trg, src, den, pot, qw);
    PROFILEEND("GPU",0);
    return;
}

template<>
void Device<GPU>::DirectStokes(const double *src, const double *den, const double *qw, 
    size_t stride, size_t n_surfs, const double *trg, size_t trg_idx_head, 
    size_t trg_idx_tail, double *pot) const
{
    std::cerr<<"The Stokes' kernel is not implemented for the double data type"<<std::endl;
    abort();
}
template<>
template<typename T>
T* Device<GPU>::ShufflePoints(T *x_in, CoordinateOrder order_in, 
    size_t stride, size_t n_surfs, T *x_out) const
{
    PROFILESTART();
    assert(x_in != x_out);
    if(order_in == AxisMajor)
        cuda_shuffle(x_in, stride, n_surfs, DIM, x_out);
    else
        cuda_shuffle(x_in, DIM, n_surfs, stride, x_out);

    PROFILEEND("GPU",0);
    return x_out;
}

template<>
template<typename T>
T Device<GPU>::Max(T *x_in, size_t length) const
{ 
    PROFILESTART();
    T m = maxGpu(x_in, length);
    PROFILEEND("GPU",0);
    return m;
}

template<>
Device<GPU>::~Device()
{
    cublasShutdown();
}
