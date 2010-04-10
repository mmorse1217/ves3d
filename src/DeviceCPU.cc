/**
 * @file   DeviceCpu.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Feb 23 15:28:14 2010
 * 
 * @brief  The implementation of the DeviceCPU class.
 */


// Device ////////////////////////////////////////////////////////////
template<typename T>
T* DeviceCPU<T>::Malloc(unsigned long int length)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::Malloc, length="<<length<<endl;
#endif
    
    return((T*) ::malloc(length * sizeof(T)));
}

template<typename T>
void DeviceCPU<T>::Free(T* ptr)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::Free"<<endl;
#endif
    
    ::free(ptr);
}

template<typename T>
T* DeviceCPU<T>::Calloc(unsigned long int num)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::Calloc, num="<<num<<endl;
#endif
    
    return((T*) ::calloc(num, sizeof(T)));
}

template<typename T>
T* DeviceCPU<T>::Memcpy(T* destination, const T* source, unsigned long int num, enum MemcpyKind kind)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::Memcpy"<<endl;
#endif
    
    if(destination == source)
    {
#ifndef NDEBUG
        cout<<"  . DeviceCPU::Memcpy, destination == source"<<endl;
#endif
        return destination;
    }
    else
    {
        return((T*) ::memcpy(destination, source, num * sizeof(T)));
    }
}

template<typename T>
T* DeviceCPU<T>::Memset(T *ptr, int value, unsigned long int num)
{
#ifndef NDEBUG
    cout<<" DeviceCPU::Memset"<<endl;
#endif
    
    return((T*) ::memset(ptr, value, num * sizeof(T)));
}

template<typename T>  
T* DeviceCPU<T>::DotProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* x_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::DotProduct"<<endl;
#endif

#ifdef PROFILING_LITE
    double ss = get_seconds();
#endif
    
    int base, resbase;
    T dot;
#pragma omp parallel for private(base, resbase, dot)
    for (int surf = 0; surf < num_surfs; surf++) {
        resbase = surf * stride;
        base = resbase * 3;
        for (int s = 0; s < stride; s++) {
            
            dot  = u_in[base + s                  ] * v_in[base + s                  ];
            dot += u_in[base + s + stride         ] * v_in[base + s + stride         ];
            dot += u_in[base + s + stride + stride] * v_in[base + s + stride + stride];
            x_out[resbase + s] = dot;
        }
    }

#ifdef PROFILING_LITE
    ss = get_seconds()-ss ;
    CpuTime::DotProduct_time +=ss;
    //cout<<"DeviceCPU::DotProduct takes (sec) : "<<ss<<endl;
#endif

    return x_out;
}

template<typename T>
T* DeviceCPU<T>::CrossProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* w_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::CrossProduct"<<endl;
#endif
    
#ifdef PROFILING
    double ss = get_seconds();
#endif

    int base, resbase, surf, s;
    T ux, uy, uz, vx, vy, vz, wx, wy, wz;
    
#pragma omp parallel for private(base, resbase, surf, s, ux, uy, uz, vx, vy, vz, wx, wy, wz)
    for (surf = 0; surf < num_surfs; surf++) {
        resbase = surf * stride;
        base = resbase * 3;
        for (s = 0; s < stride; s++) {
            ux = u_in[base + s];
            uy = u_in[base + s + stride];
            uz = u_in[base + s + stride + stride];
            vx = v_in[base + s];
            vy = v_in[base + s + stride];
            vz = v_in[base + s + stride + stride];
            wx = uy * vz - uz * vy;
            wy = uz * vx - ux * vz;
            wz = ux * vy - uy * vx;
            w_out[base + s] = wx;
            w_out[base + s + stride] = wy;
            w_out[base + s + stride + stride] = wz;
        }
    }

#ifdef PROFILING
    ss = get_seconds()-ss;
    cout<<"DeviceCPU::CrossProduct takes (sec) : "<<ss<<endl;
#endif
    return w_out;
}

template<typename T>
T* DeviceCPU<T>::Sqrt(const T* x_in, int stride, int num_surfs, T* sqrt_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::sqrt"<<endl;
#endif

#ifdef PROFILING
    double ss = get_seconds();
#endif

    int length = stride*num_surfs;
#pragma omp parallel for 
    for (int idx = 0; idx < length; idx++)
    {
        assert(x_in[idx] >= (T) 0.0);
        sqrt_out[idx] = ::sqrt(x_in[idx]);
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::Sqrt takes (sec) : "<<ss<<endl;
#endif

    return sqrt_out;
}

template<typename T>
T* DeviceCPU<T>::xInv(const T* x_in, int stride, int num_surfs, T* xInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xInv"<<endl;
#endif

#ifdef PROFILING
    double ss = get_seconds();
#endif
    
    int length = stride*num_surfs;
    T xx;    
#pragma omp parallel for private(xx)
    for (int idx = 0; idx < length; idx++)
    {
        assert(x_in[idx] != (T) 0.0);
        
        xx  = (T) 1.0;
        xx /= x_in[idx];
        xInv_out[idx] = xx;
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::xInv takes (sec) : "<<ss<<endl;
#endif

    return xInv_out;
}

template<typename T>
T* DeviceCPU<T>::xy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xy_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xy"<<endl;
#endif
    
#ifdef PROFILING_LITE
    double ss = get_seconds();
#endif
    
    int length = stride*num_surfs, idx;
    T xTy;

#pragma omp parallel for private(xTy)
    for (idx = 0; idx < length; idx++)
    {
        xTy  = x_in[idx];
        xTy *= y_in[idx];
        xy_out[idx] = xTy;
    }

#ifdef PROFILING_LITE
    ss = get_seconds()-ss ;
    CpuTime::xy_time +=ss;
    //cout<<"DeviceCPU::xy takes (sec) : "<<ss<<endl;
#endif

    return xy_out;
}

template<typename T>
T* DeviceCPU<T>::xyInv(const T* x_in, const T* y_in, int stride, int num_surfs, T* xyInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xyInv"<<endl;
#endif
    
#ifdef PROFILING
    double ss = get_seconds();
#endif
    
    int length = stride*num_surfs, idx;
    T xDy;

#pragma omp parallel for private(xDy)
    for (idx = 0; idx < length; idx++)
    {
        assert( y_in[idx] != (T) 0.0);
        xDy  = x_in[idx];
        xDy /= y_in[idx];
            
        xyInv_out[idx] = xDy;
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::xyInv takes (sec) : "<<ss<<endl;
#endif

    return xyInv_out;
}

template<typename T>
T*  DeviceCPU<T>::uyInv(const T* u_in, const T* y_in, int stride, int num_surfs, T* uyInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::uyInv"<<endl;
#endif
    
#ifdef PROFILING
    double ss = get_seconds();
#endif

    int vec, s, y_base, base, y_idx;
    T uy;

#pragma omp parallel for private(vec, s, y_base, base, y_idx, uy)
    for (vec = 0; vec < num_surfs; vec++)
    {
        y_base = vec*stride;
        base = 3*y_base;
        for (s = 0; s < stride; s++) 
        {
            y_idx = y_base + s;

            uy  = u_in[base + s];
            uy /= y_in[y_idx];
            uyInv_out[base + s] = uy;

            uy  = u_in[base + s + stride];
            uy /= y_in[y_idx];
            uyInv_out[base + s + stride] = uy;
            
            uy  = u_in[base + s + stride + stride];
            uy /= y_in[y_idx];
            uyInv_out[base + s + stride + stride] = uy;
        }
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::uyInv takes (sec) : "<<ss<<endl;
#endif

    return uyInv_out;
}

///@todo For all these methods the input and output may be the same,
///so we have temporary variables created, either make a variable or
///remove the option that in and out the same.
template<typename T>
T*  DeviceCPU<T>::axpy(T a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs , T*  axpy_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::axpy"<<endl;
#endif

#ifdef PROFILING
    double ss = get_seconds();
#endif

    T val;
    int length = stride*num_surfs;

#pragma omp parallel for private(val)
    for (int idx = 0; idx < length; idx++)
    {
        val  = a_in;
        val *= x_in[idx];
        val += y_in[idx];
        axpy_out[idx] = val;
    }
    
#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::axpy takes (sec) : "<<ss<<endl;
#endif
    
    return axpy_out;
}

template<typename T>
T*  DeviceCPU<T>::axpb(T a_in, const T*  x_in, T b_in, int stride, int num_surfs , T*  axpb_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::axpb"<<endl;
#endif

#ifdef PROFILING
    double ss = get_seconds();
#endif
    T val;
    int length = stride*num_surfs;
    
#pragma omp parallel for private(val)
    for (int idx = 0; idx < length; idx++)
    {
        val  = a_in;
        val *= x_in[idx];
        val += b_in;
        
        axpb_out[idx] = val;
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::axpb takes (sec) : "<<ss<<endl;
#endif

    return axpb_out;
}

template<typename T>
T*  DeviceCPU<T>::avpw(const T* a_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  avpw_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::avpw"<<endl;
#endif

#ifdef PROFILING
    double ss = get_seconds();
#endif
    T val;
    int base, vec, s, length = 3*stride, idx;

#pragma omp parallel for private(val, base, vec, s, idx)
    for (vec = 0; vec < num_surfs; vec++)
    {
        base = vec*length;        
        for (s = 0; s < stride; s++) {
            
            idx  = base+s;
            val  = a_in[vec];
            val *= v_in[idx];
            val += w_in[idx];
            avpw_out[idx] = val;

            idx +=stride;
            val  = a_in[vec];
            val *= v_in[idx];
            val += w_in[idx];
            avpw_out[idx] = val;
            
            idx +=stride;
            val  = a_in[vec];
            val *= v_in[idx];
            val += w_in[idx];
            avpw_out[idx] = val;
        }
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::avpw takes (sec) : "<<ss<<endl;
#endif
    
    return avpw_out;
}

template<typename T>
T*  DeviceCPU<T>::xvpw(const T* x_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  xvpw_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xvpw"<<endl;
#endif

#ifdef PROFILING
    double ss = get_seconds();
#endif
    int base, x_base, vec, s, length = 3*stride, idx, x_idx;
    T val;

#pragma omp parallel for private(base, x_base, vec, s, idx, x_idx, val)
    for (vec = 0; vec < num_surfs; vec++)
    {
        base = vec*length;
        x_base = vec*stride;
        
        for (s = 0; s < stride; s++)
        {
            idx = base+s;
            x_idx = x_base+s;
            
            val  = x_in[x_idx];
            val *= v_in[idx];
            val += w_in[idx];
            xvpw_out[idx]  = val;

            idx +=stride;
            val  = x_in[x_idx];
            val *= v_in[idx];
            val += w_in[idx];
            xvpw_out[idx]  = val;

            idx +=stride;
            val  = x_in[x_idx];
            val *= v_in[idx];
            val += w_in[idx];
            xvpw_out[idx]  = val;
        }
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::xvpw takes (sec) : "<<ss<<endl;
#endif

    return xvpw_out;
}

template<typename T>
T*  DeviceCPU<T>::xvpb(const T* x_in, const T*  v_in, T b_in, int stride, int num_surfs, T*  xvpb_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xvpb"<<endl;
#endif
    
#ifdef PROFILING_LITE
    double ss = get_seconds();
#endif

    int base, x_base, vec, s, length = 3*stride, idx, x_idx;
    T val;

#pragma omp parallel for private(val, base, x_base, vec, s, idx, x_idx)    
    for (vec = 0; vec < num_surfs; vec++)
    {
        base = vec*length;
        x_base =vec*stride;
        
        for (s = 0; s < stride; s++) {
            
            idx = base+s;
            x_idx = x_base+s;

            val  = x_in[x_idx];
            val *= v_in[idx];
            val += b_in;
            xvpb_out[idx] = val;

            idx +=stride;
            val  = x_in[x_idx];
            val *= v_in[idx];
            val += b_in;
            xvpb_out[idx] = val;

            idx +=stride;
            val  = x_in[x_idx];
            val *= v_in[idx];
            val += b_in;
            xvpb_out[idx] = val;
        }
    }

#ifdef PROFILING_LITE
    ss = get_seconds()-ss ;
    CpuTime::xvpb_time +=ss;
    //cout<<"DeviceCPU::xvpb takes (sec) : "<<ss<<endl;
#endif
    
    return xvpb_out;
}

template<typename T>
T*  DeviceCPU<T>::Reduce(const T *x_in, const T *w_in, const T *quad_w_in, int stride, int num_surfs, T  *int_x_dw)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::Reduce"<<endl;
#endif

#ifdef PROFILING
    double ss = get_seconds();
#endif

    T val, sum;
    
    if(x_in != NULL)
    {
#pragma omp parallel for private(val,sum)
        for (int ii = 0; ii < num_surfs; ++ii)
        {
            sum = 0;
            for (int jj = 0; jj < stride; ++jj) 
            {
                val  = x_in[ii * stride + jj];
                val *= w_in[ii * stride + jj];
                val *= quad_w_in[jj];
                
                sum += val;
            }
            int_x_dw[ii] = sum;
        }
    } else{
#pragma omp parallel for private(val,sum)
        for (int ii = 0; ii < num_surfs; ++ii)
        {
            sum = 0;
            for (int jj = 0; jj < stride; ++jj) 
            {
                val = w_in[ii * stride + jj];
                val *= quad_w_in[jj];
                
                sum += val;
            }

            int_x_dw[ii] = sum;
        }
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::Reduce takes (sec) : "<<ss<<endl;
#endif
    
    return int_x_dw;
}

template<typename T>
void  DeviceCPU<T>::InitializeSHT(OperatorsMats<T> &mats)
{
    assert(sht_.dft_forward == 0);
    assert(mats.fileIO_.device_ == *this);

#ifdef PROFILING
    double ss = get_seconds();
#endif

    p_ = mats.p_;
    p_up_ = mats.p_up_;

    T *dft_forward;
    T *dft_backward;
    T *dft_d1backward;
    T *dft_d2backward;

    //p version_
    dft_forward    = DeviceCPU<T>::Malloc(4 * p_ * p_); 
    dft_backward   = DeviceCPU<T>::Malloc(4 * p_ * p_); 
    dft_d1backward = DeviceCPU<T>::Malloc(4 * p_ * p_); 
    dft_d2backward = DeviceCPU<T>::Malloc(4 * p_ * p_); 
    
    sht_.InitializeBlasSht(p_, dft_forward, dft_backward, dft_d1backward, 
        dft_d2backward, mats.leg_trans_p_, mats.leg_trans_inv_p_, 
        mats.d1_leg_trans_p_, mats.d2_leg_trans_p_);

    //p_up version
    dft_forward    = DeviceCPU<T>::Malloc(4 * p_up_ * p_up_); 
    dft_backward   = DeviceCPU<T>::Malloc(4 * p_up_ * p_up_);
    dft_d1backward = DeviceCPU<T>::Malloc(4 * p_up_ * p_up_); 
    dft_d2backward = DeviceCPU<T>::Malloc(4 * p_up_ * p_up_); 

    sht_up_sample_.InitializeBlasSht(p_up_, dft_forward, dft_backward, dft_d1backward, 
        dft_d2backward, mats.leg_trans_p_up_, mats.leg_trans_inv_p_up_, 
        mats.d1_leg_trans_p_up_, mats.d2_leg_trans_p_up_);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::InitializeSHT takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
DeviceCPU<T>::DeviceCPU()
{}

template<typename T>
DeviceCPU<T>::~DeviceCPU()
{
#ifdef PROFILING_LITE
  cout<<"=========================================="<<endl;
  cout<<"=========================================="<<endl;
  cout<<"DeviceCPU::gemm  (sec) : "<<CpuTime::gemm_time<<endl;
  cout<<"DeviceCPU::stokes(sec) : "<<CpuTime::stokes_time<<endl;
  cout<<"DeviceCPU::xvpb  (sec) : "<<CpuTime::xvpb_time<<endl;
  cout<<"DeviceCPU::xy    (sec) : "<<CpuTime::xy_time<<endl;
  cout<<"DeviceCPU::Dot   (sec) : "<<CpuTime::DotProduct_time<<endl;
  cout<<"DeviceCPU::Shift (sec) : "<<CpuTime::Shift_time<<endl;
  cout<<"=========================================="<<endl;
  cout<<"=========================================="<<endl;
#endif

    Free(sht_.dft_forward); 
    Free(sht_.dft_backward); 
    Free(sht_.dft_d1backward); 
    Free(sht_.dft_d2backward); 

    // up sample
    Free(sht_up_sample_.dft_forward); 
    Free(sht_up_sample_.dft_backward); 
    Free(sht_up_sample_.dft_d1backward); 
    Free(sht_up_sample_.dft_d2backward); 
}

template<typename T>
void  DeviceCPU<T>::ShAna(const T *x_in, T *work_arr, int p, int n_funs, T *shc_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.forward(x_in, work_arr, n_funs, shc_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.forward(x_in, work_arr, n_funs, shc_out);
    }
        
        
#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShAna takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSyn(const T *shc_in, T *work_arr, int p, int n_funs, T *x_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans_inv != 0);
        sht_.backward(shc_in, work_arr, n_funs, x_out);
    }else{
        assert(sht_up_sample_.leg_trans_inv != 0);
        sht_up_sample_.backward(shc_in, work_arr, n_funs, x_out);
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::shSyn takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDu(const T *shc_in, T *work_arr, int p, int n_funs, T *xu_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.d1_leg_trans != 0);
        sht_.backward_du(shc_in, work_arr, n_funs, xu_out);
    }else{
        assert(sht_up_sample_.d1_leg_trans != 0);
        sht_up_sample_.backward_du(shc_in, work_arr, n_funs, xu_out);
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDu takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDv(const T *shc_in, T *work_arr, int p, int n_funs, T *xv_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.backward_dv(shc_in, work_arr, n_funs, xv_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.backward_dv(shc_in, work_arr, n_funs, xv_out);
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDv takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDuu(const T *shc_in, T *work_arr, int p, int n_funs, T *xuu_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.d2_leg_trans != 0);
        sht_.backward_d2u(shc_in, work_arr, n_funs, xuu_out);
    }else{
        assert(sht_up_sample_.d2_leg_trans != 0);
        sht_up_sample_.backward_d2u(shc_in, work_arr, n_funs, xuu_out);
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDuu takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDvv(const T *shc_in, T *work_arr, int p, int n_funs, T *xvv_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.leg_trans != 0);
        sht_.backward_d2v(shc_in, work_arr, n_funs, xvv_out);
    }else{
        assert(sht_up_sample_.leg_trans != 0);
        sht_up_sample_.backward_d2v(shc_in, work_arr, n_funs, xvv_out);
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDvv takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDuv(const T *shc_in, T *work_arr, int p, int n_funs, T *xuv_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif
    
    assert(p == p_ || p == p_up_);
    if ( p == p_){
        assert(sht_.d1_leg_trans != 0);
        sht_.backward_duv(shc_in, work_arr, n_funs, xuv_out);
    }else{
        assert(sht_up_sample_.d1_leg_trans != 0);
        sht_up_sample_.backward_duv(shc_in, work_arr, n_funs, xuv_out);
    }
    
#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDuv takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::AllDerivatives(const T *x_in, T *work_arr, int p, int n_funs, T *shc_x, T *Dux_out, 
    T *Dvx_out,T *Duux_out, T *Duvx_out, T *Dvvx_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

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

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::AllDerivatives takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::FirstDerivatives(const T *x_in, T *work_arr, int p, int n_funs, T *shc_x, T *Dux_out, T *Dvx_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

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

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::FirstDerivatives takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void DeviceCPU<T>::Filter(int p, int n_funs, const T *x_in, const T *alpha, T* work_arr, T *shc_out, T *x_out)
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

template<typename T>
void DeviceCPU<T>::ScaleFreqs(int p, int num_vesicles, const T * inputs, const T * alphas, T * outputs )
{
#ifndef NDEBUG
    cout<<"DeviceCPU::ScaleFreqs"<<endl;
#endif
    
    const float * inp_deb = inputs;
    float * out_deb = outputs;
    const float * alphas_deb = alphas;


    // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
    // -------- process zeroth frequency (constant) ---------
    int leg_order = p+1;
    for (int v=0; v<num_vesicles; v++)
        for (int i=0; i<leg_order; i++)
            *(outputs++) = *(inputs++) * alphas[i];
    alphas += leg_order;
    leg_order--;

    // process remaining frequencies except the last cosine
    for (; leg_order>1; leg_order--) 
    {
        // first process cosine
        for (int v=0; v<num_vesicles; v++)
            for (int i=0; i<leg_order; i++)
                *(outputs++) = *(inputs++) *alphas[i];
        alphas += leg_order;

        // then process sine
        for (int v=0; v<num_vesicles; v++)
            for (int i=0; i<leg_order; i++)
                *(outputs++) = *(inputs++) *alphas[i];
        alphas += leg_order;
    }

    // process last cosine
    for (int v=0; v<num_vesicles; v++)
        *(outputs++) = *(inputs++) * alphas[0];
    alphas += leg_order;
    leg_order--;

    assert (leg_order == 0);
    assert(inputs-inp_deb == num_vesicles*p*(p+2));
    assert(outputs-out_deb == num_vesicles*p*(p+2));
    assert(alphas-alphas_deb == p*(p+2));
}

template<typename T>
void DeviceCPU<T>::Resample(int p, int num_vesicles, int q, const T * inputs, T * outputs)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::Resample"<<endl;
#endif

    const float * inp_deb = inputs;
    float * out_deb = outputs;
    
    // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
    // -------- process zeroth frequency (constant) ---------
    int leg_order = p+1;
    int new_leg_order = q+1;
    int min_leg_order = std::min(leg_order, new_leg_order);

    for (int v=0; v<num_vesicles; v++)
    {
        for(int i=0; i<min_leg_order; i++)
            *(outputs++) = *(inputs++);
        for (int i=leg_order; i<new_leg_order; i++)
            *(outputs++) = 0;
        if (leg_order > new_leg_order)
            inputs += leg_order - new_leg_order;
    }
    leg_order--;
    new_leg_order--;
    min_leg_order--;

    // process remaining frequencies except the last cosine
    for (; min_leg_order>1; min_leg_order--,leg_order--,new_leg_order--) 
    {
        // first process cosine
        for (int v=0; v<num_vesicles; v++)
        {
            for(int i=0; i<min_leg_order; i++)
                *(outputs++) = *(inputs++);
            for (int i=leg_order; i<new_leg_order; i++)
                *(outputs++) = 0;
            if (leg_order > new_leg_order)
                inputs += leg_order - new_leg_order;
        }

        // then process sine
        for (int v=0; v<num_vesicles; v++)
        {
            for(int i=0; i<min_leg_order; i++)
                *(outputs++) = *(inputs++);
            for (int i=leg_order; i<new_leg_order; i++)
                *(outputs++) = 0;
            if (leg_order > new_leg_order)
                inputs += leg_order - new_leg_order;
        }
    }

    // process last cosine
    for (int v=0; v<num_vesicles; v++)
    {
        for(int i=0; i<min_leg_order; i++)
            *(outputs++) = *(inputs++);
        for (int i=leg_order; i<new_leg_order; i++)
            *(outputs++) = 0;
        if (leg_order > new_leg_order)
            inputs += leg_order - new_leg_order;
    }

    leg_order--;
    new_leg_order--;
    min_leg_order--;

    // assert (leg_order == 0);
 
    // if q>p all remaining coefs should be zero
    float * output_end = out_deb+num_vesicles*q*(q+2);
    assert(outputs<=output_end);

    while (outputs<output_end)
        *(outputs++) = 0;

    if (p<=q)
        assert(inputs-inp_deb == num_vesicles*p*(p+2));
    else
        assert(inputs-inp_deb < num_vesicles*p*(p+2));
    
    assert(outputs-out_deb == num_vesicles*q*(q+2));
}

template<typename T>
T* DeviceCPU<T>::gemm(const char *transA, const char *transB, const int *m, const int *n, const int *k, const T *alpha, 
    const T *A, const int *lda, const T *B, const int *ldb, const T *beta, T *C, const int *ldc)
{
    cerr<<"gemm is not implemented for this data type"<<endl;
    abort();
}



template<>
float* DeviceCPU<float>::gemm(const char *transA, const char *transB, const int *m, const int *n, const int *k, const float *alpha, 
    const float *A, const int *lda, const float *B, const int *ldb, const float *beta, float *C, const int *ldc)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::gemm"<<endl;
#endif
#ifdef PROFILING_LITE
    double ss = get_seconds();
#endif
    sgemm(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

#ifdef PROFILING_LITE
    ss = get_seconds()-ss ;
    CpuTime::gemm_time +=ss;    
//cout<<"DeviceCPU::gemm takes (sec) : "<<ss<<endl;
#endif

    return C;
}

template<typename T>
T* DeviceCPU<T>::CircShift(const T *arr_in, int n_vecs, int vec_length, int shift, T *arr_out)
{

#ifndef NDEBUG
    cout<<"DeviceCPU::CircShift"<<endl;
#endif

#ifdef PROFILING_LITE
    double ss = get_seconds();
#endif

    shift = shift%vec_length;
    if (shift<0)
        shift+=vec_length;

    int base_in, base_out;

#pragma omp parallel for private(base_in,base_out)
    for(int ii=0;ii<n_vecs; ++ii)
    {
        base_out = ii * vec_length;
        base_in  = base_out + vec_length - shift;
        for(int jj=0;jj<shift; ++jj)
            arr_out[base_out + jj] = arr_in[base_in + jj];

        base_in  = base_out;
        base_out += shift;
        for(int jj=0;jj<vec_length-shift; ++jj)
            arr_out[base_out + jj] = arr_in[base_in + jj];
    }

#ifdef PROFILING_LITE
    ss = get_seconds()-ss;
    CpuTime::Shift_time +=ss;    
    //cout<<"DeviceCPU::CircShift takes (sec) : "<<ss<<endl;
#endif

    return arr_out;
}

template<typename T>
void DeviceCPU<T>::DirectStokes(int stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const T *qw, const T *trg, const T *src, const T *den, T *pot)
{                              
#ifndef NDEBUG
    cout<<"DeviceCPU::DirectStokes"<<endl;
#endif

#ifdef PROFILING_LITE
    double ss = get_seconds();
#endif

    if(qw != NULL)
        DirectStokesKernel(stride, n_surfs, trg_idx_head, trg_idx_tail, qw, trg, src, den, pot);
    else
        DirectStokesKernel_Noqw(stride, n_surfs, trg_idx_head, trg_idx_tail, trg, src, den, pot);

#ifdef PROFILING_LITE
    ss = get_seconds()-ss;
    CpuTime::stokes_time +=ss;
    //cout<<"DeviceCPU::DirectStokes takes (sec) : "<<ss<<endl;
#endif
    return;
} 

template<>
void DeviceCPU<float>::DirectStokes(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, 
    const float *qw, const float *trg, const float *src, const float *den, float *pot)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::DirectStokes (SSE)"<<endl;
#endif

#ifdef PROFILING_LITE
    double ss = get_seconds();
#endif

    if(qw != NULL)
        DirectStokesSSE(stride, n_surfs, trg_idx_head, trg_idx_tail, qw, trg, src, den, pot);
    else
        DirectStokesKernel_Noqw(stride, n_surfs, trg_idx_head, trg_idx_tail, trg, src, den, pot);

#ifdef PROFILING_LITE
    ss = get_seconds()-ss;
    CpuTime::stokes_time += ss;
    //cout<<"DeviceCPU::DirectStokes takes (sec) : "<<ss<<endl;
#endif
    return;
}

template<typename T>
T* DeviceCPU<T>::ShufflePoints(T *x_in, CoordinateOrder order_in, int stride, int n_surfs, T *x_out)
{
    ///@todo transpose could be made in place
    assert(x_in !=x_out);

    int len = 3*stride;
    int dim1 = (order_in == AxisMajor) ? stride : 3;
    int dim2 = (order_in == AxisMajor) ? 3 : stride;

    
#pragma omp parallel for 
    for(int ss=0;ss<n_surfs;++ss)
        for(int ii=0;ii<dim1;++ii)
            for(int jj=0;jj<dim2;++jj)
                x_out[ss*len + ii*dim2+jj] = x_in[ss*len + jj*dim1+ii];
    
    return x_out;
}


template<typename T>
T DeviceCPU<T>::Max(T *x_in, int length)
{ 
    T max = *x_in;
    for(int idx = 0;idx<length;idx++)
        max = (max > x_in[idx]) ? max : x_in[idx];

    ///@todo the omp'ize version of this crashes.
    //     T *max_arr = Malloc(omp_get_num_threads());
    //     T n_threads = 0;

    // #pragma omp parallel
    //     {
    //         T max_loc = *x_in;
    // #pragma omp for
    //         for(int idx = 0;idx<length;idx++)
    //             max_loc = (max_loc > x_in[idx]) ? max_loc : x_in[idx];
    //         max_arr[omp_get_thread_num()] =  max_loc;
    //         if(omp_get_thread_num() == 0)
    //             n_threads = omp_get_num_threads();
    //     }
    
    //     T max=max_arr[0];
    //     for(int idx = 0;idx<n_threads;idx++)
    //         max = (max > max_arr[idx]) ? max : max_arr[idx];

    //    Free(max_arr);
    return(max);
}


template<typename T>
void DeviceCPU<T>::InterpSh(int p, int n_funs, const T *x_in, T* work_arr, T *shc, int q, T *x_out)
{
    ///@bug this gives segmentation fault when q<p, since x_out is small.
    ShAna(x_in, work_arr, p, n_funs, x_out);
    Resample(p, n_funs, q, x_out, shc);
    ShSyn(shc, work_arr, q, n_funs, x_out);
}
