/**
 * @file   DeviceGPU.cu
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Mon Mar  1 13:09:06 2010
 */

template <typename T>
DeviceGPU<T>::DeviceGPU(int cuda_device_id)
{
    cudaSetDevice(cuda_device_id);
}

template <typename T>
T* DeviceGPU<T>::Malloc(unsigned long int length)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::Malloc, length = "<<length<<endl;
#endif

    T* ptr = 0;
    cudaMalloc(&ptr, length * sizeof(T));
    return(ptr);
}

template <typename T>
void DeviceGPU<T>::Free(T* ptr)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::Free"<<endl;
#endif

    cudaFree(ptr);
}

template <typename T>
T* DeviceGPU<T>::Calloc(unsigned long int num)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::Calloc, num = "<<num<<endl;
#endif

    T* ptr = 0;
    cudaMalloc(&ptr, num * sizeof(T));
    cudaMemset(ptr, 0, num * sizeof(T));
    return ptr;
}

template <typename T>
T* DeviceGPU<T>::Memcpy (T* destination, const T* source, unsigned long int num, enum MemcpyKind kind)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::Memcpy, num = "<<num<<endl;
#endif

    cudaMemcpyKind cuda_kind = static_cast<cudaMemcpyKind>(kind);
    cudaMemcpy(destination, source, num * sizeof(T), cuda_kind);
    return destination;
}

template <typename T> 
T* DeviceGPU<T>::DotProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* x_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::DotProduct"<<endl;
#endif

    DotProductGpu(u_in, v_in, stride, num_surfs, x_out);
    return x_out;
}

template <typename T> 
T* DeviceGPU<T>::CrossProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* w_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::CrossProduct"<<endl;
#endif
    CrossProductGpu(u_in, v_in, stride, num_surfs, w_out); 
    return w_out;
}

template <typename T> 
T* DeviceGPU<T>::Sqrt(const T* x_in, int stride, int num_surfs, T* sqrt_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::Sqrt"<<endl;
#endif

    SqrtGpu(x_in, stride, num_surfs, sqrt_out);
    return sqrt_out;
}

template <typename T> 
T* DeviceGPU<T>::xInv(const T* x_in, int stride, int num_surfs, T* xInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::xInv"<<endl;
#endif

    InvGpu(x_in, stride, num_surfs, xInv_out);
    return xInv_out;
}

template <typename T> 
T* DeviceGPU<T>::xy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xy_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::xy"<<endl;
#endif

    xyGpu(x_in, y_in, stride, num_surfs, xy_out);
    return xy_out;
}

template <typename T> 
T* DeviceGPU<T>::xyInv(const T* x_in, const T* y_in, int stride, int num_surfs, T* xyInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::xyInv"<<endl;
#endif

    xyInvGpu(x_in, y_in, stride, num_surfs, xyInv_out);
    return xyInv_out;
}

template <typename T>
T* DeviceGPU<T>::uyInv(const T* u_in, const T* y_in, int stride, int num_surfs, T* uyInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::uyInv"<<endl;
#endif

    uyInvGpu(u_in, y_in, stride, num_surfs, uyInv_out);
    return uyInv_out;
}

template <typename T> 
T* DeviceGPU<T>::axpy(T a_in, const T* x_in, const T* y_in, int stride, int num_surfs , T* axpy_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::axpy"<<endl;
#endif

    axpyGpu(a_in, x_in, y_in, stride, num_surfs, axpy_out);
    return axpy_out;
}

template <typename T> 
T* DeviceGPU<T>::axpb(T a_in, const T*  x_in, T b_in, int stride, int num_surfs , T*  axpb_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::axpb"<<endl;
#endif

    axpbGpu(a_in, x_in, b_in, stride, num_surfs, axpb_out);
    return axpb_out;
}

template <typename T> 
T* DeviceGPU<T>::avpw(const T* a_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  avpw_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::avpw"<<endl;
#endif

    cout<<"GPU AVPW"<<endl;
    return avpw_out;
}

template <typename T> 
T* DeviceGPU<T>::xvpw(const T* x_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  xvpw_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::xvpw"<<endl;
#endif

    xvpwGpu(v_in, x_in, w_in, stride, num_surfs, xvpw_out);
    return xvpw_out;
}

template <typename T> 
T* DeviceGPU<T>::xvpb(const T* x_in, const T*  v_in, T b_in, int stride, int num_surfs, T*  xvpb_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::xvpb"<<endl;
#endif

    xvpbGpu(v_in, x_in, b_in, stride, num_surfs, xvpb_out);
    return xvpb_out;
}

template <typename T> 
T* DeviceGPU<T>::Reduce(const T *x_in, const T *w_in, const T *quad_w_in, int stride, int num_surfs, T  *int_x_dw)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::Reduce"<<endl;
#endif

    ReduceGpu(x_in, w_in, quad_w_in, stride, num_surfs, int_x_dw);
    return int_x_dw;
}

template <typename T> 
T* DeviceGPU<T>::gemm(const char *transA, const char *transB, const int *m, const int *n, const int *k, const T *alpha, 
    const T *A, const int *lda, const T *B, const int *ldb, const T *beta, T *C, const int *ldc)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::gemm"<<endl;
#endif

    ///@bug transA may not be consistent
    cublasSgemm(*transA, *transB, *m, *n, *k, *alpha, A, *lda, B, *ldb, *beta, C, *ldc); 
    return C;
}

template <typename T> 
T* DeviceGPU<T>::CircShift(const T *arr_in, int n_vecs, int vec_length, int shift, T *arr_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::CircShift"<<endl;
#endif

    CircShiftGpu(arr_in, n_vecs, vec_length, shift, arr_out);
    return arr_out;
}


template <typename T>
void DeviceGPU<T>::DirectStokes(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, 
    const T *qw, const T *trg, const T *src, const T *den, T *pot)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::DirectStokes"<<endl;
#endif
    cuda_stokes(stride, n_surfs, trg_idx_head, trg_idx_tail, trg, src, den, pot, qw);
}

template <typename T> 
T* DeviceGPU<T>::ShufflePoints(T *x_in, CoordinateOrder order_in, int stride, int n_surfs, T *x_out)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::ShufflePoints"<<endl;
#endif

    Memcpy(x_out, x_in, 3 * stride * n_surfs, MemcpyDeviceToDevice);
    int dim = 3;
    if(order_in == AxisMajor)
        cuda_shuffle(x_out, stride, n_surfs, dim);
    else
        cuda_shuffle(x_out, dim, n_surfs, stride);

    return x_out;
}


template <typename T>
T DeviceGPU<T>::Max(T *x_in, int length)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::Max"<<endl;
#endif

    cout<<"GPU MAX"<<endl;
    return(-1);
}

template <typename T>
void DeviceGPU<T>::InitializeSHT(OperatorsMats<T> &mats)
{
#ifndef NDEBUG
    cout<<"DeviceGPU::InitializeSHT"<<endl;
#endif

    assert(sht_.dft_forward == 0);
    assert(mats.fileIO_.device_ == *this);
    p_ = mats.p_;
    p_up_ = mats.p_up_;

    T *dft_forward;
    T *dft_backward;
    T *dft_d1backward;
    T *dft_d2backward;

    //p version_
    dft_forward    = Malloc(4 * p_ * p_); 
    dft_backward   = Malloc(4 * p_ * p_); 
    dft_d1backward = Malloc(4 * p_ * p_); 
    dft_d2backward = Malloc(4 * p_ * p_); 

    sht_.InitializeCudaSht(p_, dft_forward, dft_backward, dft_d1backward, 
        dft_d2backward, mats.leg_trans_p_, mats.leg_trans_inv_p_, 
        mats.d1_leg_trans_p_, mats.d2_leg_trans_p_);

    //p_up version
    dft_forward    = Malloc(4 * p_up_ * p_up_); 
    dft_backward   = Malloc(4 * p_up_ * p_up_); 
    dft_d1backward = Malloc(4 * p_up_ * p_up_); 
    dft_d2backward = Malloc(4 * p_up_ * p_up_); 

    sht_up_sample_.InitializeCudaSht(p_up_, dft_forward, dft_backward, dft_d1backward, 
        dft_d2backward, mats.leg_trans_p_up_, mats.leg_trans_inv_p_up_, 
        mats.d1_leg_trans_p_up_, mats.d2_leg_trans_p_up_);
    
    cublasInit();
}

template <typename T>
DeviceGPU<T>::~DeviceGPU()
{
    Free(sht_.dft_forward); 
    Free(sht_.dft_backward); 
    Free(sht_.dft_d1backward); 
    Free(sht_.dft_d2backward); 
 
    // up sample
    Free(sht_up_sample_.dft_forward); 
    Free(sht_up_sample_.dft_backward); 
    Free(sht_up_sample_.dft_d1backward); 
    Free(sht_up_sample_.dft_d2backward); 
    
    // for gemm 
    cublasShutdown();
}

template <typename T>
void DeviceGPU<T>::ShAna(const T *x_in, T *work_arr, int p, int n_funs, T *shc_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.forward(x_in, work_arr, n_funs, shc_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.forward(x_in, work_arr, n_funs, shc_out);
    }
}

template <typename T>
void DeviceGPU<T>::ShSyn(const T *shc_in, T *work_arr, int p, int n_funs, T *y_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans_inv != 0);
        sht_.backward(shc_in, work_arr, n_funs, y_out);
    }else{
        assert(sht_up_sample_.leg_trans_inv != 0);
        sht_up_sample_.backward(shc_in, work_arr, n_funs, y_out);
    }
}

template <typename T>
void DeviceGPU<T>::ShSynDu(const T *shc_in, T *work_arr, int p, int n_funs, T *xu_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.d1_leg_trans != 0);
        sht_.backward_du(shc_in, work_arr, n_funs, xu_out);
    }else{
        assert(sht_up_sample_.d1_leg_trans != 0);
        sht_up_sample_.backward_du(shc_in, work_arr, n_funs, xu_out);
    }
}

template <typename T>
void DeviceGPU<T>::ShSynDv(const T *shc_in, T *work_arr, int p, int n_funs, T *xv_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.backward_dv(shc_in, work_arr, n_funs, xv_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.backward_dv(shc_in, work_arr, n_funs, xv_out);
    }
}

template <typename T>
void DeviceGPU<T>::ShSynDuu(const T *shc_in, T *work_arr, int p, int n_funs, T *xuu_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.d2_leg_trans != 0);
        sht_.backward_d2u(shc_in, work_arr, n_funs, xuu_out);
    }else{
        assert(sht_up_sample_.d2_leg_trans != 0);
        sht_up_sample_.backward_d2u(shc_in, work_arr, n_funs, xuu_out);
    }
}

template <typename T>
void DeviceGPU<T>::ShSynDvv(const T *shc_in, T *work_arr, int p, int n_funs, T *xvv_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.backward_d2v(shc_in, work_arr, n_funs, xvv_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.backward_d2v(shc_in, work_arr, n_funs, xvv_out);
    }
}

template <typename T>
void DeviceGPU<T>::ShSynDuv(const T *shc_in, T *work_arr, int p, int n_funs, T *xuv_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.d1_leg_trans != 0);
        sht_.backward_duv(shc_in, work_arr, n_funs, xuv_out);
    }else{
        assert(sht_up_sample_.d1_leg_trans != 0);
        sht_up_sample_.backward_duv(shc_in, work_arr, n_funs, xuv_out);
    }
}

template <typename T>
void DeviceGPU<T>::AllDerivatives(const T *x_in, T *work_arr, int p, int n_funs, T* shc_x, T *Dux_out, T *Dvx_out, 
    T *Duux_out, T *Duvx_out, T *Dvvx_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.forward(     x_in , work_arr, n_funs, shc_x);
        sht_.backward_du( shc_x, work_arr, n_funs, Dux_out);
        sht_.backward_dv( shc_x, work_arr, n_funs, Dvx_out);
        sht_.backward_d2u(shc_x, work_arr, n_funs, Duux_out);
        sht_.backward_d2v(shc_x, work_arr, n_funs, Dvvx_out);
        sht_.backward_duv(shc_x, work_arr, n_funs, Duvx_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.forward(     x_in , work_arr, n_funs, shc_x);
        sht_up_sample_.backward_du( shc_x, work_arr, n_funs, Dux_out);
        sht_up_sample_.backward_dv( shc_x, work_arr, n_funs, Dvx_out);
        sht_up_sample_.backward_d2u(shc_x, work_arr, n_funs, Duux_out);
        sht_up_sample_.backward_d2v(shc_x, work_arr, n_funs, Dvvx_out);
        sht_up_sample_.backward_duv(shc_x, work_arr, n_funs, Duvx_out);
    }
}

template <typename T>
void DeviceGPU<T>::FirstDerivatives(const T *x_in, T *work_arr, int p, int n_funs, T* shc_x, T *Dux_out, T *Dvx_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.forward(     x_in , work_arr, n_funs, shc_x);
        sht_.backward_du( shc_x, work_arr, n_funs, Dux_out);
        sht_.backward_dv( shc_x, work_arr, n_funs, Dvx_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.forward(     x_in , work_arr, n_funs, shc_x);
        sht_up_sample_.backward_du( shc_x, work_arr, n_funs, Dux_out);
        sht_up_sample_.backward_dv( shc_x, work_arr, n_funs, Dvx_out);
    }
}

template <typename T>
void DeviceGPU<T>::Filter(int p, int n_funs, const T *x_in, const T *alpha, T* work_arr, T *shc_out, T *x_out)
{
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.forward(x_in, work_arr, n_funs, shc_out);
        ScaleFreqs(p, n_funs, shc_out, alpha, shc_out);
        sht_.backward(shc_out, work_arr, n_funs, x_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.forward(x_in, work_arr, n_funs, shc_out);
        ScaleFreqs(p, n_funs, shc_out, alpha, shc_out);
        sht_up_sample_.backward(shc_out, work_arr, n_funs, x_out);
    }
}

template <typename T>
void DeviceGPU<T>::ScaleFreqs(int p, int n_funs, const T *inputs, const T *alphas, T *outputs)
{
    ScaleFreqsGpu(p, n_funs, inputs, alphas, outputs);
}

template <typename T>
void DeviceGPU<T>::Resample(int p, int n_funs, int q, const T *shc_p, T *shc_q)
{
    ResampleGpu(p, n_funs, q, shc_p, shc_q);
}


