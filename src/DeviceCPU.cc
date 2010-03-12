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

    int ii, jj, idx = 0;
    T val;

#pragma omp parallel for private(val,ii, jj, idx)
    for (ii = 0; ii < num_surfs; ++ii)
    {
        int_x_dw[ii] = 0;
        for (jj = 0; jj < stride; ++jj) 
        {
            val  = x_in[idx];
            val *= w_in[idx++];
            val *= quad_w_in[jj];

            int_x_dw[ii] += val;
        }
    }

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::Reduce takes (sec) : "<<ss<<endl;
#endif
    
    return int_x_dw;
}

template<typename T>
void  DeviceCPU<T>::InitializeSHT(int p, char *leg_trans_fname,
    char *leg_trans_inv_fname, char *d1_leg_trans_fname, 
    char *d2_leg_trans_fname)
{
    assert(sht_.dft_forward == 0);
    
#ifdef PROFILING
    double ss = get_seconds();
#endif

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

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::InitializeSHT takes (sec) : "<<ss<<endl;
#endif
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
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(sht_.leg_trans != 0);
    sht_.forward(x_in, work_arr, n_funs, shc_out);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShAna takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSyn(const T *shc_in, T *work_arr, int n_funs, T *x_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(sht_.leg_trans_inv != 0);
    sht_.backward(shc_in, work_arr, n_funs, x_out);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::shSyn takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDu(const T *shc_in, T *work_arr, int n_funs, T *xu_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(sht_.d1_leg_trans != 0);
    sht_.backward_du(shc_in, work_arr, n_funs, xu_out);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDu takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDv(const T *shc_in, T *work_arr, int n_funs, T *xv_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(sht_.leg_trans != 0);
    sht_.backward_dv(shc_in, work_arr, n_funs, xv_out);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDv takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDuu(const T *shc_in, T *work_arr, int n_funs, T *xuu_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(sht_.d2_leg_trans != 0);
    sht_.backward_d2u(shc_in, work_arr, n_funs, xuu_out);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDuu takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::ShSynDvv(const T *shc_in, T *work_arr, int n_funs, T *xvv_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif
    assert(sht_.leg_trans != 0);
    sht_.backward_d2v(shc_in, work_arr, n_funs, xvv_out);
}

template<typename T>
void  DeviceCPU<T>::ShSynDuv(const T *shc_in, T *work_arr, int n_funs, T *xuv_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(sht_.d1_leg_trans != 0);
    sht_.backward_duv(shc_in, work_arr, n_funs, xuv_out);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::ShSynDuv takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::AllDerivatives(const T *x_in, T *work_arr, int n_funs, T *shc_x, T *Dux_out, 
    T *Dvx_out,T *Duux_out, T *Duvx_out, T *Dvvx_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(sht_.leg_trans != 0);
    sht_.forward(     x_in , work_arr, n_funs, shc_x);
    sht_.backward_du( shc_x, work_arr, n_funs, Dux_out);
    sht_.backward_dv( shc_x, work_arr, n_funs, Dvx_out);
    sht_.backward_d2u(shc_x, work_arr, n_funs, Duux_out);
    sht_.backward_d2v(shc_x, work_arr, n_funs, Dvvx_out);
    sht_.backward_duv(shc_x, work_arr, n_funs, Duvx_out);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::AllDerivatives takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void  DeviceCPU<T>::FirstDerivatives(const T *x_in, T *work_arr, int n_funs, T *shc_x, T *Dux_out, T *Dvx_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif

    assert(sht_.leg_trans != 0);
    sht_.forward(x_in, work_arr, n_funs, shc_x);
    sht_.backward_du(shc_x, work_arr, n_funs, Dux_out);
    sht_.backward_dv(shc_x, work_arr, n_funs, Dvx_out);

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"DeviceCPU::FirstDerivatives takes (sec) : "<<ss<<endl;
#endif
}

template<typename T>
void DeviceCPU<T>::Filter(int p, int n_funs, const T *x_in, const T *alpha, T* work_arr, T *shc_out, T *x_out)
{
    assert(sht_.leg_trans != 0);

    sht_.forward(x_in, work_arr, n_funs, shc_out);
    ScaleFreqs(p, n_funs, shc_out, alpha, shc_out);
    sht_.backward(shc_out, work_arr, n_funs, x_out);
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

