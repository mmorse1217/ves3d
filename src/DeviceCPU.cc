/**
 * @file   DeviceCpu.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Feb 23 15:28:14 2010
 * 
 * @brief  The implementation of the DeviceCPU class.
 */

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

    int base, resbase;
    T dot;

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
    return x_out;
}

template<typename T>
T* DeviceCPU<T>::CrossProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* w_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::CrossProduct"<<endl;
#endif
    
    int base, resbase, surf, s;
    T ux, uy, uz, vx, vy, vz, wx, wy, wz;

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
    return w_out;
}

template<typename T>
T* DeviceCPU<T>::Sqrt(const T* x_in, int stride, int num_surfs, T* sqrt_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::sqrt"<<endl;
#endif

    int length = stride*num_surfs;
    for (int idx = 0; idx < length; idx++)
    {
        assert(x_in[idx] >= (T) 0.0);
        sqrt_out[idx] = ::sqrt(x_in[idx]);
    }
    return sqrt_out;
}

template<typename T>
T* DeviceCPU<T>::xInv(const T* x_in, int stride, int num_surfs, T* xInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xInv"<<endl;
#endif

    int length = stride*num_surfs;
    for (int idx = 0; idx < length; idx++)
    {
        assert(x_in[idx] != (T) 0.0);
        xInv_out[idx] = 1.0/x_in[idx];
    }

    return xInv_out;
}

template<typename T>
T* DeviceCPU<T>::xy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xy_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xy"<<endl;
#endif
    
    int length = stride*num_surfs, idx;

    for (idx = 0; idx < length; idx++)
        xy_out[idx] = x_in[idx] * y_in[idx];

    return xy_out;
}

template<typename T>
T* DeviceCPU<T>::xyInv(const T* x_in, const T* y_in, int stride, int num_surfs, T* xyInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xyInv"<<endl;
#endif
    
    int length = stride*num_surfs, idx;

    for (idx = 0; idx < length; idx++)
    {
        assert( y_in[idx] != (T) 0.0);
        xyInv_out[idx] = x_in[idx]/y_in[idx];
    }
    return xyInv_out;
}

template<typename T>
T*  DeviceCPU<T>::uyInv(const T* u_in, const T* y_in, int stride, int num_surfs, T* uyInv_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::uyInv"<<endl;
#endif
    
    int vec, s, y_base, base, y_idx;

    for (vec = 0; vec < num_surfs; vec++)
    {
        y_base = vec*stride;
        base = 3*y_base;
        for (s = 0; s < stride; s++) 
        {
            y_idx = y_base + s;
            uyInv_out[base                  + s ]  = u_in[base                   + s ] / y_in[y_idx];
            uyInv_out[base+ stride          + s ]  = u_in[base + stride          + s ] / y_in[y_idx];
            uyInv_out[base+ stride + stride + s ]  = u_in[base + stride + stride + s ] / y_in[y_idx];
        }
    }
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

    T val;
    int length = stride*num_surfs;
    for (int idx = 0; idx < length; idx++)
        axpy_out[idx] = a_in*x_in[idx] + y_in[idx];
    
    return axpy_out;
}

template<typename T>
T*  DeviceCPU<T>::axpb(T a_in, const T*  x_in, T b_in, int stride, int num_surfs , T*  axpb_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::axpb"<<endl;
#endif

    int length = stride*num_surfs;
    for (int idx = 0; idx < length; idx++)
        axpb_out[idx] = a_in*x_in[idx] + b_in;

    return axpb_out;
}

template<typename T>
T*  DeviceCPU<T>::xvpw(const T* x_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  xvpw_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xvpw"<<endl;
#endif

    int base, x_base, vec, s, length = 3*stride, idx, x_idx;

    for (vec = 0; vec < num_surfs; vec++)
    {
        base = vec*length;
        x_base =vec*stride;
        
        for (s = 0; s < stride; s++) {
            
            idx = base+s;
            x_idx = x_base+s;

            xvpw_out[idx]  = x_in[x_idx] * v_in[idx] + w_in[idx];

            idx +=stride;
            xvpw_out[idx]  = x_in[x_idx] * v_in[idx] + w_in[idx];

            idx +=stride;
            xvpw_out[idx]  = x_in[x_idx] * v_in[idx] + w_in[idx];
        }
    }
    return xvpw_out;
}

template<typename T>
T*  DeviceCPU<T>::xvpb(const T* x_in, const T*  v_in, T b_in, int stride, int num_surfs, T*  xvpb_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::xvpb"<<endl;
#endif
    int base, x_base, vec, s, length = 3*stride, idx, x_idx;

    for (vec = 0; vec < num_surfs; vec++)
    {
        base = vec*length;
        x_base =vec*stride;
        
        for (s = 0; s < stride; s++) {
            
            idx = base+s;
            x_idx = x_base+s;

            xvpb_out[idx]  = x_in[x_idx] * v_in[idx];
            xvpb_out[idx] += b_in;

            idx +=stride;
            xvpb_out[idx]  = x_in[x_idx] * v_in[idx];
            xvpb_out[idx] += b_in;

            idx +=stride;
            xvpb_out[idx]  = x_in[x_idx] * v_in[idx];
            xvpb_out[idx] += b_in;
        }
    }

    return xvpb_out;
}

template<typename T>
void  DeviceCPU<T>::InitializeSHT(int p, char *leg_trans_fname,
    char *leg_trans_inv_fname, char *d1_leg_trans_fname, 
    char *d2_leg_trans_fname)
{
    assert(dft_forward == 0);

     T *dft_forward;
    T *dft_backward;
    T *dft_d1backward;
    T *dft_d2backward;
    T *leg_trans;
    T *leg_trans_inv;
    T *d1_leg_trans;
    T *d2_leg_trans;

    dft_forward    = DeviceCPU<T>::Malloc(4 * p * p); 
    dft_backward   = DeviceCPU<T>::Malloc(4 * p * p); 
    dft_d1backward = DeviceCPU<T>::Malloc(4 * p * p); 
    dft_d2backward = DeviceCPU<T>::Malloc(4 * p * p); 

    leg_trans      = DeviceCPU<T>::Malloc((p + 1) * (p+1) * (p +2)); 
    leg_trans_inv  = DeviceCPU<T>::Malloc((p + 1) * (p+1) * (p +2)); 
    d1_leg_trans   = DeviceCPU<T>::Malloc((p + 1) * (p+1) * (p +2)); 
    d2_leg_trans   = DeviceCPU<T>::Malloc((p + 1) * (p+1) * (p +2)); 

    sht_.InitializeBlasSht(p, leg_trans_fname, leg_trans_inv_fname, 
        d1_leg_trans_fname, d2_leg_trans_fname, 
        dft_forward, dft_backward, dft_d1backward, dft_d2backward, 
        leg_trans, leg_trans_inv, d1_leg_trans, d2_leg_trans);
}

template<typename T>
DeviceCPU<T>::DeviceCPU(){}

template<typename T>
DeviceCPU<T>::~DeviceCPU()
{
    if(sht_.dft_forward != 0)
        Free(sht_.dft_forward); 
    
    if(sht_.dft_backward != 0)
        Free(sht_.dft_backward); 

    if(sht_.dft_d1backward != 0)
        Free(sht_.dft_d1backward); 
    
    if(sht_.dft_d2backward != 0)
        Free(sht_.dft_d2backward); 
    
    if(sht_.leg_trans != 0)
        Free(sht_.leg_trans); 
    
    if(sht_.leg_trans_inv != 0)
        Free(sht_.leg_trans_inv); 
    
    if(sht_.d1_leg_trans != 0)
        Free(sht_.d1_leg_trans); 

    if(sht_.d2_leg_trans != 0)
        Free(sht_.d2_leg_trans); 
}

template<typename T>
void  DeviceCPU<T>::ShAna(const T *x_in, T *work_arr, int n_funs, T *shc_out)
{
    assert(sht_.leg_trans != 0);
    sht_.forward(x_in, work_arr, n_funs, shc_out);
}

template<typename T>
void  DeviceCPU<T>::ShSyn(const T *shc_in, T *work_arr, int n_funs, T *x_out)
{
    assert(sht_.leg_trans_inv != 0);
    sht_.backward(shc_in, work_arr, n_funs, x_out);
}

template<typename T>
void  DeviceCPU<T>::ShSynDu(const T *shc_in, T *work_arr, int n_funs, T *xu_out)
{
    assert(sht_.d1_leg_trans != 0);
    sht_.backward_du(shc_in, work_arr, n_funs, xu_out);
}

template<typename T>
void  DeviceCPU<T>::ShSynDv(const T *shc_in, T *work_arr, int n_funs, T *xv_out)
{
    assert(sht_.leg_trans != 0);
    sht_.backward_dv(shc_in, work_arr, n_funs, xv_out);
}

template<typename T>
void  DeviceCPU<T>::ShSynDuu(const T *shc_in, T *work_arr, int n_funs, T *xuu_out)
{
    assert(sht_.d2_leg_trans != 0);
    sht_.backward_d2u(shc_in, work_arr, n_funs, xuu_out);
}

template<typename T>
void  DeviceCPU<T>::ShSynDvv(const T *shc_in, T *work_arr, int n_funs, T *xvv_out)
{
    assert(sht_.leg_trans != 0);
    sht_.backward_d2v(shc_in, work_arr, n_funs, xvv_out);
}

template<typename T>
void  DeviceCPU<T>::ShSynDuv(const T *shc_in, T *work_arr, int n_funs, T *xuv_out)
{
    assert(sht_.d1_leg_trans != 0);
    sht_.backward_duv(shc_in, work_arr, n_funs, xuv_out);
}

template<typename T>
void  DeviceCPU<T>::AllDerivatives(const T *x_in, T *work_arr, int n_funs, T *shc_x, T *Dux_out, 
    T *Dvx_out,T *Duux_out, T *Duvx_out, T *Dvvx_out)
{
    assert(sht_.leg_trans != 0);
    sht_.forward(     x_in , work_arr, n_funs, shc_x);
    sht_.backward_du( shc_x, work_arr, n_funs, Dux_out);
    sht_.backward_dv( shc_x, work_arr, n_funs, Dvx_out);
    sht_.backward_d2u(shc_x, work_arr, n_funs, Duux_out);
    sht_.backward_d2v(shc_x, work_arr, n_funs, Dvvx_out);
    sht_.backward_duv(shc_x, work_arr, n_funs, Duvx_out);
}

template<typename T>
void  DeviceCPU<T>::FirstDerivatives(const T *x_in, T *work_arr, int n_funs, T *shc_x, T *Dux_out, T *Dvx_out)
{
    assert(sht_.leg_trans != 0);
    sht_.forward(x_in, work_arr, n_funs, shc_x);
    sht_.backward_du(shc_x, work_arr, n_funs, Dux_out);
    sht_.backward_dv(shc_x, work_arr, n_funs, Dvx_out);
}

