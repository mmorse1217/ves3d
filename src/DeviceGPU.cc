template<>
Device<GPU>::Device(int device_id, Error_t *err)
{
    if(err != NULL)
        switch ( cudaSetDevice(device_id) )
        {
            case cudaSuccess:
                *err = Success;
                break;
            case cudaErrorInvalidDevice:
                *err = InvalidDevice;
                break;
            case cudaErrorSetOnActiveProcess:
                *err = SetOnActiveDevice;
                break;
        }
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
    PROFILEEND("GPU", 5 * n_vecs * stride);
    return x_out;
}

template<>
template<typename T>  
T* Device<GPU>::CrossProduct(const T* u_in, const T* v_in, 
    size_t stride, size_t num_surfs, T* w_out) const
{
    PROFILESTART();
    assert(DIM==3);
    CrossProductGpu(u_in, v_in, stride, num_surfs, w_out); 
    PROFILEEND("GPU", 9 * stride * num_surfs);
    return w_out;
}

template<>
template<typename T>
T* Device<GPU>::Sqrt(const T* x_in, size_t length, T* sqrt_out) const
{
    PROFILESTART();
    SqrtGpu(x_in, length, sqrt_out);
    PROFILEEND("GPU", length);
    return sqrt_out;
}


template<>
template<typename T>
T* Device<GPU>::ax(const T* a, const T* x, size_t length, 
    size_t n_vecs, T* ax_out)  const
{   
    PROFILESTART();
    axGpu(a, x, length, n_vecs, ax_out);
    PROFILEEND("GPU", length * n_vecs);
    return(ax_out);
}

template<>
template<typename T>
T* Device<GPU>::xy(const T* x_in, const T* y_in, size_t length, T* xy_out) const
{
    PROFILESTART();
    xyGpu(x_in, y_in, length, xy_out);
    PROFILEEND("GPU", length);
    return xy_out;
}

template<>
template<typename T>
T* Device<GPU>::xyInv(const T* x_in, const T* y_in, 
    size_t length, T* xyInv_out) const
{
    PROFILESTART();
    if(x_in==NULL)
        InvGpu(y_in, length, xyInv_out);
    else
        xyInvGpu(x_in, y_in, length, xyInv_out);

    PROFILEEND("GPU",length);
    return xyInv_out;
}

template<>
template<typename T>
T*  Device<GPU>::uyInv(const T* u_in, const T* y_in, 
    size_t stride, size_t num_surfs, T* uyInv_out) const
{
    PROFILESTART();
    assert(u_in!=NULL);
    uyInvGpu(u_in, y_in, stride, num_surfs, uyInv_out);
    PROFILEEND("GPU", stride * num_surfs);
    return uyInv_out;
}

template<>
template<typename T>
T*  Device<GPU>::axpy(T a_in, const T*  x_in, const T*  y_in, 
    size_t length, T*  axpy_out) const
{
    PROFILESTART();
    assert(x_in != NULL);
    if(y_in !=NULL)
        axpyGpu(a_in, x_in, y_in, length, axpy_out);
    else
        axpbGpu(a_in, x_in, (T) 0.0, length, axpy_out);
        
    PROFILEEND("GPU", 2 * length);
    
    return axpy_out;
}

template<>
template<typename T>
T* Device<GPU>::apx(const T* a_in, const T* x_in, size_t stride, 
        size_t n_subs, T* apx_out) const
{
    PROFILESTART();
    apxGpu(a_in, x_in, stride, n_subs, apx_out);
    PROFILEEND("GPU", n_subs * stride);

    return(apx_out);
}

template<>
template<typename T>
T*  Device<GPU>::avpw(const T* a_in, const T*  v_in, const T*  w_in, 
    size_t stride, size_t num_surfs, T*  avpw_out) const
{
    PROFILESTART();
    if(w_in !=NULL)
    {
        avpwGpu(a_in, v_in, w_in, stride, num_surfs, avpw_out);
        PROFILEEND("GPU", 6 * num_surfs * stride);
    }
    else
    {
        avpwGpu(a_in, v_in, stride, num_surfs, avpw_out);
        PROFILEEND("GPU", 3 * num_surfs * stride);
    }
    return avpw_out;
}

template<>
template<typename T>
T*  Device<GPU>::xvpw(const T* x_in, const T*  v_in, const T*  w_in, 
    size_t stride, size_t num_surfs, T*  xvpw_out) const
{
    PROFILESTART();
    if(w_in !=NULL)
    {
        xvpwGpu(v_in, x_in, w_in, stride, num_surfs, xvpw_out);
        PROFILEEND("GPU", 6 * stride * num_surfs);
    }
    else
    {
        xvpbGpu(v_in, x_in, (T) 0.0, stride, num_surfs, xvpw_out);
        PROFILEEND("GPU", 6 * stride * num_surfs);
    }
    return xvpw_out;
}

template<>
template<typename T>
T*  Device<GPU>::Reduce(const T *x_in, const int x_dim, const T *w_in, 
    const T *quad_w_in, size_t stride, size_t num_surfs, T  *int_x_dw) const
{
    PROFILESTART();
    ReduceGpu(x_in, x_dim, w_in, quad_w_in, stride, num_surfs, int_x_dw);
    PROFILEEND("GPU", ((x_in == NULL) ? 2 : 3 * x_dim) * num_surfs * stride);
    return int_x_dw;
}

template<>
template<typename T>
T* Device<GPU>::gemm(const char *transA, const char *transB, 
    const int *m, const int *n, const int *k, const T *alpha, 
    const T *A, const int *lda, const T *B, const int *ldb, 
    const T *beta, T *C, const int *ldc) const
{
    PROFILESTART();
    cugemm(transA, transB, m, n, k, alpha, A, lda, B, 
        ldb, beta, C, ldc); 
    cudaThreadSynchronize();
    PROFILEEND("GPU",(double) 2* (*k) * (*n) * (*m) + *(beta) * (*n) * (*m));
    return C;
}

template<>
template<typename T>
void Device<GPU>::DirectStokes(const T *src, const T *den,
    const T *qw, size_t stride, size_t n_surfs, const T *trg,
    size_t trg_idx_head, size_t trg_idx_tail, T *pot) const
{
    PROFILESTART();
    cuda_stokes(stride, n_surfs, trg_idx_head, trg_idx_tail, 
        trg, src, den, pot, qw);
    PROFILEEND("GPU",((qw == NULL) ? 32 : 35) * n_surfs * stride * (trg_idx_tail - trg_idx_head));
    return;
}

template<>
template<typename T>
T Device<GPU>::MaxAbs(T *x_in, size_t length) const
{ 
    PROFILESTART();
    T m = maxGpu(x_in, length);
    PROFILEEND("GPU", 0);
    return m;
}

template<>
template<typename T>
T* Device<GPU>::Transpose(const T *in, size_t height, 
    size_t width, T *out) const
{
    PROFILESTART();
    cu_trans(out, in, width, height);
    PROFILEEND("GPU", 0);
    return(out);
}

template<>
template<typename T>
T Device<GPU>::AlgebraicDot(const T* x, const T* y, size_t length) const
{
    PROFILESTART();
    T dot = AlgebraicDotGpu(x, y, length);
    PROFILEEND("GPU", length);
    return(dot);
}

template<>
Device<GPU>::~Device()
{
    cublasShutdown();
}

template<>
template<typename T>
void Device<GPU>::AggregateRotation(int sh_order, int n_vec,
    const int* n_sub, const T* mat, const T** vec, T** res) const
{
    PROFILESTART();
    
    int nlat = gridDimOf(sh_order).first;
    int nlong = gridDimOf(sh_order).second;
    int np = nlat * nlong;
    
    T alpha(1), beta(0);
        
    cudaStream_t *stream = new cudaStream_t[nlong];
    T** rot_mat = new T*[nlong];
    
    for (int jj = 0; jj < nlong; ++jj) 
    {
        cudaStreamCreate(&stream[jj]);
        rot_mat[jj] = (T*) this->Malloc(np * np * sizeof(T));
    }

    for (int jj=0; jj< nlong; ++jj)
    {
        PermuteGpu(mat, sh_order, jj, rot_mat[jj], stream[jj]);
        cublasSetKernelStream(stream[jj]);
        for(int ii=0; ii<n_vec; ++ii)
        {
            int nsub(n_sub[ii]);
            this->gemm("N", "N", &np, &nsub, &np, &alpha, rot_mat[jj], &np, 
                vec[ii], &np, &beta, res[n_vec * jj + ii], &np);
        }
    }
    cudaThreadSynchronize();
    
    for (int jj = 0; jj < nlong; ++jj) 
        this->Free(rot_mat[jj]);
    
    delete[] stream;
    delete[] rot_mat;

    PROFILEEND("",0);
}
