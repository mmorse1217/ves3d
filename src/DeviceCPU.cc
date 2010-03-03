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
T*  DeviceCPU<T>::axpy(T a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs , T*  axpy_out)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::axpy"<<endl;
#endif

    int length = stride*num_surfs;
    for (int idx = 0; idx < length; idx++)
    {
        axpy_out[idx]  = a_in*x_in[idx];
        axpy_out[idx] += y_in[idx];
    }
    
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
    {
        axpb_out[idx] = a_in*x_in[idx];
        axpb_out[idx] += b_in;
    }

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

            xvpw_out[idx]  = x_in[x_idx] * v_in[idx];
            xvpw_out[idx] += w_in[idx];

            idx +=stride;
            xvpw_out[idx]  = x_in[x_idx] * v_in[idx];
            xvpw_out[idx] += w_in[idx];

            idx +=stride;
            xvpw_out[idx]  = x_in[x_idx] * v_in[idx];
            xvpw_out[idx] += w_in[idx];
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
