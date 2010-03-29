/**
 * @file   DeviceCpu.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Feb 23 15:28:14 2010
 * 
 * @brief  The implementation of the DeviceCPU class.
 */

// OperatorsMats //////////////////////////////////////////////////////////////////////

template <typename T>
OperatorsMats<T>::OperatorsMats(int p_in, int p_up_in, bool readFromFile) : 
    p_(p_in), p_up_(p_up_in)
{
    int np = 2 * p_ * ( p_ + 1);
    struct DeviceCPU<T> dev; ///@bug This should be moved outside AND IO class should be fixed
    DataIO<T> myIO(dev,"",0); 

    quad_weights_ = (T*) malloc(np * sizeof(T));
    all_rot_mats_ = (T*) malloc( np * np * (p_ + 1) * sizeof(T));
    sing_quad_weights_ = (T*) malloc(np * sizeof(T));
    w_sph_ = (T*) malloc(np * sizeof(T));

    int leg_size = (p_ + 1) * (p_+1) * (p_ +2);
    leg_trans_p_     = (T*) malloc( leg_size * sizeof(T)); 
    leg_trans_inv_p_ = (T*) malloc( leg_size * sizeof(T));
    d1_leg_trans_p_  = (T*) malloc( leg_size * sizeof(T));
    d2_leg_trans_p_  = (T*) malloc( leg_size * sizeof(T));
     
    //p_up_
    quad_weights_p_up_ = (T*) malloc(2 * p_up_ *(p_up_ + 1) * sizeof(T));
    
    int leg_size_up = (p_up_ + 1) * (p_up_+1) * (p_up_ +2);	
    leg_trans_p_up_     = (T*) malloc( leg_size_up * sizeof(T)); 
    leg_trans_inv_p_up_ = (T*) malloc( leg_size_up * sizeof(T));
    d1_leg_trans_p_up_  = (T*) malloc( leg_size_up * sizeof(T));
    d2_leg_trans_p_up_  = (T*) malloc( leg_size_up * sizeof(T));
    
    char fname[300];
    if(readFromFile)
    {
        sprintf(fname,"%s/precomputed/quad_weights_%u_single.txt",getenv("VES3D_DIR"),p_);
        //sprintf(fname,"precomputed/quad_weights_%u_single.txt",p_);
        myIO.ReadData(fname, np, quad_weights_);

        sprintf(fname,"%s/precomputed/sing_quad_weights_%u_single.txt",getenv("VES3D_DIR"),p_);
        //sprintf(fname,"precomputed/sing_quad_weights_%u_single.txt",p_);
        myIO.ReadData(fname, np, sing_quad_weights_);
    
        sprintf(fname,"%s/precomputed/all_rot_mats_%u_single.txt",getenv("VES3D_DIR"),p_);
        //sprintf(fname,"precomputed/all_rot_mats_%u_single.txt",p_);
        myIO.ReadData(fname, np * np *(p_ + 1), all_rot_mats_);
    
        sprintf(fname,"%s/precomputed/w_sph_%u_single.txt",getenv("VES3D_DIR"),p_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_);
        myIO.ReadData(fname, np, w_sph_);

        sprintf(fname,"%s/precomputed/legTrans%u_single.txt",getenv("VES3D_DIR"),p_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_);
        myIO.ReadData(fname, leg_size, leg_trans_p_);

        sprintf(fname,"%s/precomputed/legTransInv%u_single.txt",getenv("VES3D_DIR"),p_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_);
        myIO.ReadData(fname, leg_size, leg_trans_inv_p_);

        sprintf(fname,"%s/precomputed/d1legTrans%u_single.txt",getenv("VES3D_DIR"),p_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_);
        myIO.ReadData(fname, leg_size, d1_leg_trans_p_);

        sprintf(fname,"%s/precomputed/d2legTrans%u_single.txt",getenv("VES3D_DIR"),p_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_);
        myIO.ReadData(fname, leg_size, d2_leg_trans_p_);

        //p_up_
        sprintf(fname,"%s/precomputed/quad_weights_%u_single.txt",getenv("VES3D_DIR"),p_up_);
        //sprintf(fname,"precomputed/quad_weights_%u_single.txt",p_up_);
        myIO.ReadData(fname, 2 * p_up_ *(p_up_ + 1), quad_weights_p_up_);
        
        sprintf(fname,"%s/precomputed/legTrans%u_single.txt",getenv("VES3D_DIR"),p_up_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_up_);
        myIO.ReadData(fname, leg_size, leg_trans_p_up_);

        sprintf(fname,"%s/precomputed/legTransInv%u_single.txt",getenv("VES3D_DIR"),p_up_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_up_);
        myIO.ReadData(fname, leg_size, leg_trans_inv_p_up_);

        sprintf(fname,"%s/precomputed/d1legTrans%u_single.txt",getenv("VES3D_DIR"),p_up_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_up_);
        myIO.ReadData(fname, leg_size, d1_leg_trans_p_up_);

        sprintf(fname,"%s/precomputed/d2legTrans%u_single.txt",getenv("VES3D_DIR"),p_up_);
        //sprintf(fname,"precomputed/w_sph_%u_single.txt",p_up_);
        myIO.ReadData(fname, leg_size, d2_leg_trans_p_up_);
    }
}

template <typename T>
OperatorsMats<T>::~OperatorsMats()
{
    free(quad_weights_);
    free(all_rot_mats_);
    free(sing_quad_weights_);
    free(w_sph_);

    free(leg_trans_p_);
    free(leg_trans_inv_p_);
    free(d1_leg_trans_p_);
    free(d2_leg_trans_p_);
    
    free(quad_weights_p_up_);
    free(leg_trans_p_up_);
    free(leg_trans_inv_p_up_);
    free(d1_leg_trans_p_up_);
    free(d2_leg_trans_p_up_);

}

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
T* DeviceCPU<T>::DotProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* x_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::DotProduct"<<endl;
#endif

#ifdef PROFILING
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

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::DotProduct takes (sec) : "<<ss<<endl;
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
    
#ifdef PROFILING
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

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::xy takes (sec) : "<<ss<<endl;
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
    
#ifdef PROFILING
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

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::xvpb takes (sec) : "<<ss<<endl;
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
void DeviceCPU<T>::ScaleFreqs(int p, int n_funs, const T *shc_in, const T *alpha, T *shc_out)
{
    const T *inp_deb = shc_in;
    const T *alpha_deb = alpha;
    T *out_deb = shc_out;
    
    // we have even-order real dft; this means we don't have first sine
    // (sine of zero frequency) and last sine (sine of half-order
    // frequency) -------- process zeroth frequency (constant) ---------
    int leg_order = p+1;
    for (int v=0; v<n_funs; v++)
        for (int i=0; i<leg_order; i++)
            *(shc_out++) = *(shc_in++) * alpha[i];
    alpha += leg_order;
    leg_order--;

    // process remaining frequencies except the last cosine
    for (; leg_order>1; leg_order--) 
    {
        // first process cosine
        for (int v=0; v<n_funs; v++)
            for (int i=0; i<leg_order; i++)
                *(shc_out++) = *(shc_in++) *alpha[i];
        alpha += leg_order;
        
        // then process sine
        for (int v=0; v<n_funs; v++)
            for (int i=0; i<leg_order; i++)
                *(shc_out++) = *(shc_in++) *alpha[i];
        alpha += leg_order;
    }
    
    // process last cosine
    for (int v=0; v<n_funs; v++)
        for (int i=0; i<leg_order; i++)
            *(shc_out++) = *(shc_in++) * alpha[i];
    alpha += leg_order;
    leg_order--;
    
    assert(leg_order == 0);
    assert(shc_in-inp_deb == n_funs*p*(p+2));
    assert(shc_out-out_deb == n_funs*p*(p+2));
    assert(alpha-alpha_deb == p*(p+2));
}

template<typename T>
void DeviceCPU<T>::Resample(int p, int n_funs, int q, const T *shc_p, T *shc_q)
{
    const T * inp_deb = shc_p;
    T * out_deb = shc_q;

    // we have even-order real dft; this means we don't have first sine
    // (sine of zero frequency) and last sine (sine of half-order
    // frequency) -------- process zeroth frequency (constant) ---------
    int leg_order = p+1;
    int new_leg_order = q+1;
    int min_leg_order = (leg_order < new_leg_order) ? leg_order : new_leg_order;

    for (int v=0; v<n_funs; v++)
    {
        for(int i=0; i<min_leg_order; i++)
            *(shc_q++) = *(shc_p++);
        for (int i=leg_order; i<new_leg_order; i++)
            *(shc_q++) = 0;
        if (leg_order > new_leg_order)
            shc_p += leg_order - new_leg_order;
    }
    leg_order--;
    new_leg_order--;
    min_leg_order--;

    // process remaining frequencies except the last cosine
    for (; min_leg_order>1; min_leg_order--,leg_order--,new_leg_order--) 
    {
        // first process cosine
        for (int v=0; v<n_funs; v++)
        {
            for(int i=0; i<min_leg_order; i++)
                *(shc_q++) = *(shc_p++);
            for (int i=leg_order; i<new_leg_order; i++)
                *(shc_q++) = 0;
            if (leg_order > new_leg_order)
                shc_p += leg_order - new_leg_order;
        }

        // then process sine
        for (int v=0; v<n_funs; v++)
        {
            for(int i=0; i<min_leg_order; i++)
                *(shc_q++) = *(shc_p++);
            for (int i=leg_order; i<new_leg_order; i++)
                *(shc_q++) = 0;
            if (leg_order > new_leg_order)
                shc_p += leg_order - new_leg_order;
        }
    }

    // process last cosine
    for (int v=0; v<n_funs; v++)
        *(shc_q++) = *(shc_p++);

    leg_order--;
    new_leg_order--;
    min_leg_order--;

    assert (min_leg_order == 0);
 
    // if q>p all remaining coefs should be zero
    T * output_end = out_deb+n_funs*q*(q+2);
    assert(shc_q<=output_end);

    while (shc_q<output_end)
        *(shc_q++) = 0;

    if (p<=q)
        assert(shc_p-inp_deb == n_funs*p*(p+2));
    else
        assert(shc_p-inp_deb < n_funs*p*(p+2));
    
    assert(shc_q-out_deb == n_funs*q*(q+2));
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
#ifdef PROFILING
    double ss = get_seconds();
#endif
    sgemm(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::gemm takes (sec) : "<<ss<<endl;
#endif
    return C;
}

template<typename T>
T* DeviceCPU<T>::CircShift(const T *arr_in, int n_vecs, int vec_length, int shift, T *arr_out)
{

#ifndef NDEBUG
    cout<<"DeviceCPU::CircShift"<<endl;
#endif

#ifdef PROFILING
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

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::CircShift takes (sec) : "<<ss<<endl;
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

#ifdef PROFILING
    double ss = get_seconds();
#endif

    T tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;
    
#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc)
    for (int vt=0; vt<n_surfs; vt++)
    {
        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            px = 0;
            py = 0;
            pz = 0;
            
            tx=trg[3*vt*stride +                   trg_idx];
            ty=trg[3*vt*stride + stride +          trg_idx];
            tz=trg[3*vt*stride + stride + stride + trg_idx];
            
            for (int s=0; s<stride; s++)
            {
                dx=src[3*stride*vt +                   s]-tx;
                dy=src[3*stride*vt + stride +          s]-ty;
                dz=src[3*stride*vt + stride + stride + s]-tz;

                invR = dx*dx;
                invR+= dy*dy;
                invR+= dz*dz;
                
                if (invR!=0)
                    invR = 1.0/sqrt(invR);
            
                cpx = den[3*stride*vt +                   s] * qw[s]; 
                cpy = den[3*stride*vt + stride +          s] * qw[s]; 
                cpz = den[3*stride*vt + stride + stride + s] * qw[s]; 
                
                cc  = dx*cpx;
                cc += dy*cpy;
                cc += dz*cpz;
                cc *= invR;
                cc *= invR;

                cpx += cc*dx;
                cpy += cc*dy;
                cpz += cc*dz;
                
                px += cpx*invR;
                py += cpy*invR;
                pz += cpz*invR;
            }
            pot[3*vt*stride +                  trg_idx] = px;
            pot[3*vt*stride + stride +         trg_idx] = py;
            pot[3*vt*stride + stride +stride + trg_idx] = pz;
        }
    }

#ifdef PROFILING
    ss = get_seconds()-ss;
    cout<<"DeviceCPU::DirectStokes takes (sec) : "<<ss<<endl;
#endif
    return;
} 

#define IDEAL_ALIGNMENT 16
#define SIMD_LEN (IDEAL_ALIGNMENT / sizeof(float))

// template<>
// void DeviceCPU<float>::DirectStokes(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, 
//     const float *qw, const float *trg, const float *src, const float *den, float *pot)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::DirectStokes (SSE, float)"<<endl;
// #endif

// #ifdef PROFILING
//     double ss = get_seconds();
// #endif

//     //#ifdef __SSE2__ 
//     ///@todo check for availability of SSE

//     if (stride%4) // necessary for proper alignment of sources
//         abort();
    
//     float aux_arr[3*SIMD_LEN+3]; 
//     float *tempvalx; 
//     float *tempvaly; 
//     float *tempvalz; 
//     size_t residual = size_t(aux_arr)%IDEAL_ALIGNMENT;
//     if (residual)  // if aux_arr is misaligned
//         tempvalx = aux_arr + (IDEAL_ALIGNMENT - residual);
//     else tempvalx = aux_arr;
//     if (size_t(tempvalx)%IDEAL_ALIGNMENT)  // for debugging
//         abort();
//     tempvaly=tempvalx+SIMD_LEN;
//     tempvalz=tempvaly+SIMD_LEN;
  
//     //#pragma omp parallel for private(aux_arr)
//     ///@todo add the openmp instructions
//     for (int vt=0; vt<n_surfs; vt++)
//     {
//         float p[3]={0,0,0};
//         float tx=trg[3*vt*stride            + trg_idx];
//         float ty=trg[3*vt*stride +   stride + trg_idx];
//         float tz=trg[3*vt*stride + 2*stride + trg_idx];

//         residual = size_t(src+ 3*stride*vt)%IDEAL_ALIGNMENT;
//         if (residual)
//             residual = IDEAL_ALIGNMENT - residual;
        
//         // Handle start data if it is not 16-byte aligned
//         size_t s;
//         // residual = stride;
//         for (s=0; s<residual; s++)
//         {
//             float dX_reg=src[3*stride*vt+           s]-tx;
//             float dY_reg=src[3*stride*vt+  stride + s]-ty;
//             float dZ_reg=src[3*stride*vt+2*stride + s]-tz;
            
//             float invR = (dX_reg*dX_reg+dY_reg*dY_reg+dZ_reg*dZ_reg);
//             if (invR!=0)
//                 invR = 1.0/sqrt(invR);
            
//             float cur_pot_x = den[3*stride*vt +           s] * qw[s];
//             float cur_pot_y = den[3*stride*vt +  stride + s] * qw[s];
//             float cur_pot_z = den[3*stride*vt +2*stride + s] * qw[s];

//             float tmp_scalar = (dX_reg*cur_pot_x + dY_reg*cur_pot_y + dZ_reg*cur_pot_z)*invR*invR;
//             cur_pot_x += tmp_scalar*dX_reg;
//             cur_pot_y += tmp_scalar*dY_reg;
//             cur_pot_z += tmp_scalar*dZ_reg;
            
//             p[0] += cur_pot_x*invR;
//             p[1] += cur_pot_y*invR;
//             p[2] += cur_pot_z*invR;
//         }

//         __m128 txi = _mm_load1_ps (&tx);
//         __m128 tyi = _mm_load1_ps (&ty);
//         __m128 tzi = _mm_load1_ps (&tz);

//         // for (int s=0; s<stride; s++)
//         __m128 tempx;
//         __m128 tempy;
//         __m128 tempz;

//         tempx = _mm_setzero_ps();
//         tempy = _mm_setzero_ps();
//         tempz = _mm_setzero_ps();
        
//         // Load and calculate in groups of SIMD_LEN
//         size_t loop_limit = stride-SIMD_LEN;
//         for (; s <= loop_limit; s += SIMD_LEN) {
//             __m128 sxj = _mm_load_ps (src+3*stride*vt+         s);
//             __m128 syj = _mm_load_ps (src+3*stride*vt+  stride+s);
//             __m128 szj = _mm_load_ps (src+3*stride*vt+2*stride+s);
            
//             // this could be vectorized assuming den and q are 16-byte aligned
//             __m128 sdenx = _mm_set_ps (
//                 den[3*stride*vt +s+3] * qw[s+3],
//                 den[3*stride*vt +s+2] * qw[s+2],
//                 den[3*stride*vt +s+1] * qw[s+1],
//                 den[3*stride*vt +s] * qw[s]);

//             __m128 sdeny = _mm_set_ps (
//                 den[3*stride*vt+stride +s+3] * qw[s+3],
//                 den[3*stride*vt+stride +s+2] * qw[s+2],
//                 den[3*stride*vt+stride +s+1] * qw[s+1],
//                 den[3*stride*vt+stride +s] * qw[s]);

//             __m128 sdenz = _mm_set_ps (
//                 den[3*stride*vt+2*stride +s+3] * qw[s+3],
//                 den[3*stride*vt+2*stride +s+2] * qw[s+2],
//                 den[3*stride*vt+2*stride +s+1] * qw[s+1],
//                 den[3*stride*vt+2*stride +s] * qw[s]
//                                        );
            
//             //       __m128 sdenx = _mm_load_ps (src+3*stride*vt+         s);
//             //       __m128 sdeny = _mm_load_ps (src+3*stride*vt+  stride+s);
//             //       __m128 sdenz = _mm_load_ps (src+3*stride*vt+2*stride+s);
            
//             __m128 dX, dY, dZ;
//             __m128 dR2;
//             __m128 S;

//             dX = _mm_sub_ps(txi , sxj);
//             dY = _mm_sub_ps(tyi , syj);
//             dZ = _mm_sub_ps(tzi , szj);

//             sxj = _mm_mul_ps(dX, dX); 
//             syj = _mm_mul_ps(dY, dY);
//             szj = _mm_mul_ps(dZ, dZ);

//             dR2 = _mm_add_ps(sxj, syj);
//             dR2 = _mm_add_ps(szj, dR2);

//             __m128 zero = _mm_setzero_ps ();
//             __m128 is_zero = _mm_cmpeq_ps(dR2, zero);

//             // S = _mm_rsqrt_ps(dR2);
//             const __m128 approx = _mm_rsqrt_ps( dR2 );
//             const __m128 muls = _mm_mul_ps(_mm_mul_ps(dR2, approx), approx);
//             const __m128 three = _mm_set1_ps (3.0f);
//             const __m128 half4 = _mm_set1_ps (0.5f);
//             S = _mm_mul_ps(_mm_mul_ps(half4, approx), _mm_sub_ps(three, muls) );
//             S = _mm_andnot_ps (is_zero, S);
            
//             __m128 dotx = _mm_mul_ps (dX, sdenx);
//             __m128 doty = _mm_mul_ps (dY, sdeny);
//             __m128 dotz = _mm_mul_ps (dZ, sdenz);

//             __m128 dot_sum = _mm_add_ps (dotx, doty);
//             dot_sum = _mm_add_ps (dot_sum, dotz);

//             dot_sum = _mm_mul_ps (dot_sum, S);
//             dot_sum = _mm_mul_ps (dot_sum, S);

//             dotx = _mm_mul_ps (dot_sum, dX);
//             doty = _mm_mul_ps (dot_sum, dY);
//             dotz = _mm_mul_ps (dot_sum, dZ);

//             sdenx = _mm_add_ps (sdenx, dotx);
//             sdeny = _mm_add_ps (sdeny, doty);
//             sdenz = _mm_add_ps (sdenz, dotz);

//             sdenx = _mm_mul_ps (sdenx, S);
//             sdeny = _mm_mul_ps (sdeny, S);
//             sdenz = _mm_mul_ps (sdenz, S);

//             tempx = _mm_add_ps (sdenx, tempx);
//             tempy = _mm_add_ps (sdeny, tempy);
//             tempz = _mm_add_ps (sdenz, tempz);

//         }
        
//         _mm_store_ps(tempvalx, tempx); 
//         _mm_store_ps(tempvaly, tempy); 
//         _mm_store_ps(tempvalz, tempz); 
        
//         for (size_t k = 0; k < SIMD_LEN; k++) {
//             p[0] += tempvalx[k];
//             p[1] += tempvaly[k];
//             p[2] += tempvalz[k];
//         }
        
//         if (s!=size_t(stride))
//             abort();
        
//         pot[3*vt*stride +            trg_idx] = p[0];
//         pot[3*vt*stride +   stride + trg_idx] = p[1];
//         pot[3*vt*stride + 2*stride + trg_idx] = p[2];
//     }

// #ifdef PROFILING
//     ss = get_seconds()-ss ;
//     cout<<"DeviceCPU::DirectStokes (SSE - float) takes (sec) : "<<ss<<endl;
// #endif
//     return;
// }

template<typename T>
T* DeviceCPU<T>::ShufflePoints(T *x_in, CoordinateOrder order_in, int stride, int n_surfs, T *x_out)
{

    assert(x_out != x_in);
        
    int idx_in, idx_out;
    T x, y, z;
    if(order_in == AxisMajor)
    {
#pragma omp parallel for private(idx_in, idx_out, x, y, z)
        for(int ii=0;ii<n_surfs;++ii)
        {
            idx_in  = 3*ii*stride;
            idx_out = idx_in-1;
                
            for(int jj=0;jj<stride;++jj)
            {
                x = x_in[idx_in                   ];
                y = x_in[idx_in + stride          ];
                z = x_in[idx_in + stride + stride ];
                idx_in++;
                    
                x_out[++idx_out] = x;
                x_out[++idx_out] = y;
                x_out[++idx_out] = z;
            }
        }
    }
    else
    {
#pragma omp parallel for private(idx_in, idx_out, x, y, z)
        for(int ii=0;ii<n_surfs;++ii)
        {
            idx_out = 3*ii*stride;
            idx_in  = idx_out-1;
                
            for(int jj=0;jj<stride;++jj)
            {
                x = x_in[++idx_in];
                y = x_in[++idx_in];
                z = x_in[++idx_in];
                    
                x_out[idx_out                   ];
                x_out[idx_out + stride          ];
                x_out[idx_out + stride + stride ];
                idx_out++;
            }
        }
    }
    
    return x_out;
}
