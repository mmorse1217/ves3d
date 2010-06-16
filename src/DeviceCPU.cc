/**
 * @file   DeviceCPU.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Feb 23 15:28:14 2010
 * 
 * @brief  The implementation of the Device class.
 */

inline pair<int, int> gridDimOf(int sh_order)
{
    return(make_pair(sh_order + 1, 2 * sh_order));
}

std::ostream& operator<<(std::ostream& output, const enum DeviceType &DT)
{
    switch (DT)
    {
        case CPU:
            output<<"CPU";
            break;
        case GPU:
            output<<"GPU";
            break;
    }
    
    return output;
}    

std::ostream& operator<<(std::ostream& output, const enum MemcpyKind &MK)
{
    switch (MK)
    {
        case MemcpyHostToHost:
            output<<"MemcpyHostToHost";
            break;

        case MemcpyHostToDevice:
            output<<"MemcpyHostToDevice";
            break;

        case MemcpyDeviceToHost:
            output<<"MemcpyDeviceToHost";
            break;

        case MemcpyDeviceToDevice:
            output<<"MemcpyDeviceToDevice";
            break;
    }
    
    return output;
}    

std::ostream& operator<<(std::ostream& output, const enum DeviceError &err)
{
    switch (err)
    {
        case Success:
            output<<"DeviceError: Success";
            break;
        case InvalidDevice:
            output<<"DeviceError: InvalidDevice";
            break;
        case SetOnActiveDevice:
            output<<"DeviceError: SetOnActiveDevice";
            break;
    }
    
    return output;
}    

template<enum DeviceType DT>
Device<DT>::Device(int device_id, enum DeviceError *err)
{
    if(err!=0) *err = Success;

}

template<enum DeviceType DT>
Device<DT>::~Device()
{}

template<enum DeviceType DT>
void* Device<DT>::Malloc(size_t length) const
{
    PROFILESTART();
    void* ptr = ::malloc(length);
    PROFILEEND("CPU",0);
    return(ptr);
}

template<enum DeviceType DT>
void Device<DT>::Free(void* ptr) const
{
    PROFILESTART();
    ::free(ptr);
    ptr = 0;
    PROFILEEND("CPU",0);
}

template<enum DeviceType DT>
void* Device<DT>::Calloc(size_t num, size_t size) const
{
    PROFILESTART();
    void * ptr = ::calloc(num, size);
    PROFILEEND("CPU",0);
    return(ptr);
}

template<enum DeviceType DT>
void* Device<DT>::Memcpy(void* destination, const void* source, 
    size_t num, enum MemcpyKind kind) const
{
    PROFILESTART();
    ::memcpy(destination, source, num);
    PROFILEEND("CPU",0);
    return(destination);
}

template<enum DeviceType DT>
void* Device<DT>::Memset(void *ptr, int value, size_t num) const
{
    PROFILESTART();
    ::memset(ptr, value, num);
    PROFILEEND("CPU",0);
    return(ptr);
}

template<enum DeviceType DT>
template<typename T>  
T* Device<DT>::DotProduct(const T* u_in, const T* v_in, size_t stride, 
    size_t n_vecs, T* x_out) const
{
    PROFILESTART();
    int base, resbase;
    T dot;
#pragma omp parallel for private(base, resbase, dot)
    for (int vv = 0; vv < n_vecs; vv++) {
        resbase = vv * stride;
        base = resbase * DIM;
        for (int s = 0; s < stride; s++) {
            dot = 0.0;
            for(int dd=0; dd<DIM;++dd)
                dot  += u_in[base + s + dd*stride] * v_in[base + s + dd*stride];
            
            x_out[resbase + s] = dot;
        }
    }

    PROFILEEND("CPU",0);
    return x_out;
}

template<enum DeviceType DT>
template<typename T>  
T* Device<DT>::CrossProduct(const T* u_in, const T* v_in, size_t stride, size_t num_surfs, T* w_out) const
{
    PROFILESTART();
    assert(DIM==3);
    
#pragma omp parallel 
    {
        T u[DIM], v[DIM], w[DIM];
        size_t base, resbase, surf, s;
        
#pragma omp for
        for (surf = 0; surf < num_surfs; surf++) {
            resbase = surf * stride;
            base = resbase * DIM;
            for(s = 0; s < stride; s++) {
                for(int dd=0;dd<DIM;++dd)
                {
                    u[dd] = u_in[base + s + dd * stride];
                    v[dd] = v_in[base + s + dd * stride];
                }

                w[0] = u[1] * v[2] - u[2] * v[1];
                w[1] = u[2] * v[0] - u[0] * v[2];
                w[2] = u[0] * v[1] - u[1] * v[0];
                
                for(int dd=0;dd<DIM;++dd)
                    w_out[base + s + dd * stride] = w[dd];
            }
        }
    }

    PROFILEEND("CPU",0);
    return w_out;
}

template<enum DeviceType DT>
template<typename T>
T* Device<DT>::Sqrt(const T* x_in, size_t length, T* sqrt_out) const
{
    PROFILESTART();
#pragma omp parallel for 
    for (int idx = 0; idx < length; idx++)
    {
        assert(x_in[idx] >= (T) 0.0);
        sqrt_out[idx] = ::sqrt(x_in[idx]);
    }
    PROFILEEND("CPU",0);
    return sqrt_out;
}

template<enum DeviceType DT>
template<typename T>
T* Device<DT>::ax(const T* a, const T* x, size_t stride, size_t n_vecs, T* ax_out) const
{
    PROFILESTART();
        
#pragma omp parallel for
    for(size_t n = 0;  n < n_vecs; ++n)
        for (size_t idx = 0; idx < stride; ++idx)
            ax_out[n*stride + idx]  = a[idx] * x[n * stride + idx];
    

    PROFILEEND("CPU",0);
    return ax_out;
}

template<enum DeviceType DT>
template<typename T>
T* Device<DT>::xy(const T* x_in, const T* y_in, size_t length, T* xy_out) const
{
    PROFILESTART();
    T xTy;

#pragma omp parallel for private(xTy)
    for (size_t idx = 0; idx < length; idx++)
    {
        xTy  = x_in[idx];
        xTy *= y_in[idx];
        xy_out[idx] = xTy;
    }

    PROFILEEND("CPU",0);
    return xy_out;
}

template<enum DeviceType DT>
template<typename T>
T* Device<DT>::xyInv(const T* x_in, const T* y_in, size_t length, T* xyInv_out) const
{
    PROFILESTART();
    T xDy;

    if(x_in == NULL)
    {
#pragma omp parallel for private(xDy)
        for (size_t idx = 0; idx < length; idx++)
        {
            assert( y_in[idx] != (T) 0.0);
            xDy  = 1.0;
            xDy /= y_in[idx];
            
            xyInv_out[idx] = xDy;
        }
    }
    else
    {
#pragma omp parallel for private(xDy)
        for (size_t idx = 0; idx < length; idx++)
        {
            assert( y_in[idx] != (T) 0.0);
            xDy  = x_in[idx];
            xDy /= y_in[idx];
            
            xyInv_out[idx] = xDy;
        }
    }

    PROFILEEND("CPU",0);
    return xyInv_out;
}

template<enum DeviceType DT>
template<typename T>
T*  Device<DT>::uyInv(const T* u_in, const T* y_in, size_t stride, size_t num_surfs, T* uyInv_out) const
{
    PROFILESTART();
    assert(u_in!=NULL);
    
    size_t y_base, base, y_idx;
    T uy;
    
#pragma omp parallel for private(y_base, base, y_idx, uy)
    for (size_t vec = 0; vec < num_surfs; vec++)
    {
        y_base = vec*stride;
        base = DIM*y_base;
        for (size_t s = 0; s < stride; s++) 
        {
            y_idx = y_base + s;
            for(int dd=0;dd<DIM;++dd)
            {
                uy  = u_in[base + s + dd*stride];
                uy /= y_in[y_idx];
                uyInv_out[base + s + dd*stride] = uy;
            }
        }
    }
    
    PROFILEEND("CPU",0);
    return uyInv_out;
}

template<enum DeviceType DT>
template<typename T>
T*  Device<DT>::axpy(T a_in, const T*  x_in, const T*  y_in, size_t length, T*  axpy_out) const
{
    PROFILESTART();
    assert(x_in != NULL);
    
    T val;
    
    if(y_in !=NULL)
    {
#pragma omp parallel for private(val)
        for (size_t idx = 0; idx < length; idx++)
        {
            val  = a_in;
            val *= x_in[idx];
            val += y_in[idx];
            axpy_out[idx] = val;
        }
    }
    else
    {
#pragma omp parallel for private(val)
        for (size_t idx = 0; idx < length; idx++)
        {
            val  = a_in;
            val *= x_in[idx];
            axpy_out[idx] = val;
        }
    }    

    PROFILEEND("CPU",0);
    return axpy_out;
}

template<enum DeviceType DT>
template<typename T>
T*  Device<DT>::apx(T* a_in, const T* x_in, size_t stride, 
    size_t n_subs, T* apx_out) const
{
    PROFILESTART();
    assert(a_in != NULL);
    assert(x_in != NULL);
    
#pragma omp parallel for
    for (size_t ii = 0; ii < n_subs; ++ii)
        for(size_t jj=0; jj < stride; ++jj)
            apx_out[ii * stride + jj] = a_in[ii] + x_in[ii * stride + jj];
           

    PROFILEEND("CPU",0);
    return apx_out;
}

template<enum DeviceType DT>
template<typename T>
T*  Device<DT>::avpw(const T* a_in, const T*  v_in, const T*  w_in, size_t stride, size_t num_surfs, T*  avpw_out) const
{
    PROFILESTART();
    T val;
    size_t base, vec, s, length = DIM*stride, idx;

    if(w_in !=NULL)
    {
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
    }
    else
    {
#pragma omp parallel for private(val, base, vec, s, idx)
        for (vec = 0; vec < num_surfs; vec++)
        {
            base = vec*length;        
            for (s = 0; s < stride; s++) {
            
                idx  = base+s;
                val  = a_in[vec];
                val *= v_in[idx];
                avpw_out[idx] = val;

                idx +=stride;
                val  = a_in[vec];
                val *= v_in[idx];
                avpw_out[idx] = val;
            
                idx +=stride;
                val  = a_in[vec];
                val *= v_in[idx];
                avpw_out[idx] = val;
            }
        }
    }

    PROFILEEND("CPU",0);
    return avpw_out;
}

template<enum DeviceType DT>
template<typename T>
T*  Device<DT>::xvpw(const T* x_in, const T*  v_in, const T*  w_in, size_t stride, size_t num_surfs, T*  xvpw_out) const
{
    PROFILESTART();
    size_t base, x_base, vec, s, length = DIM*stride, idx, x_idx;
    T val;

    if(w_in !=NULL)
    {
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

    }
    else
    {

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
            xvpw_out[idx]  = val;

            idx +=stride;
            val  = x_in[x_idx];
            val *= v_in[idx];
            xvpw_out[idx]  = val;

            idx +=stride;
            val  = x_in[x_idx];
            val *= v_in[idx];
            xvpw_out[idx]  = val;
        }
    }
    }

    PROFILEEND("CPU",0);
    return xvpw_out;
}

template<enum DeviceType DT>
template<typename T>
T*  Device<DT>::Reduce(const T *x_in, const int x_dim, const T *w_in, const T *quad_w_in, 
    const size_t stride, const size_t ns, T *x_dw) const
{
    PROFILESTART();
    T val, sum;
    
    if(x_in != NULL)
    {
#pragma omp parallel for private(val,sum)
        for (size_t ii = 0; ii < ns; ++ii)
        {
            int wbase = ii * stride;
            for(int kk = 0; kk < x_dim; ++kk)
            {
                sum = 0;
                int xbase = ii * x_dim * stride + kk * stride;
                
                for (size_t jj = 0; jj < stride; ++jj) 
                {
                    val  = x_in[xbase + jj];
                    val *= w_in[wbase + jj];
                    val *= quad_w_in[jj];
                    
                    sum += val;
                }
                x_dw[ii*x_dim + kk] = sum;
            }
        }
    } else{
#pragma omp parallel for private(val,sum)
        for (size_t ii = 0; ii < ns; ++ii)
        {
            sum = 0;
            for (size_t jj = 0; jj < stride; ++jj) 
            {
                val = w_in[ii * stride + jj];
                val *= quad_w_in[jj];
                
                sum += val;
            }

           x_dw[ii] = sum;
        }
    }

    PROFILEEND("CPU",0);
    return x_dw;
}

template<enum DeviceType DT>
float* Device<DT>::gemm(const char *transA, const char *transB, 
    const int *m, const int *n, const int *k, const float *alpha, 
    const float *A, const int *lda, const float *B, const int *ldb, 
    const float *beta, float *C, const int *ldc) const
{
    PROFILESTART();
    sgemm(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    PROFILEEND("CPUs", (double) 2* (*k) * (*n) * (*m) + *(beta) * (*n) * (*m));
    return C;
}

template<enum DeviceType DT>
double* Device<DT>::gemm(const char *transA, const char *transB, 
    const int *m, const int *n, const int *k, const double *alpha, 
    const double *A, const int *lda, const double *B, const int *ldb, 
    const double *beta, double *C, const int *ldc) const
{
    PROFILESTART();
    dgemm(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    PROFILEEND("CPUd",(double) 2* (*k) * (*n) * (*m) + *(beta) * (*n) * (*m));
    return C;
}

template<enum DeviceType DT>
void Device<DT>::DirectStokes(const double *src, const double *den, const double *qw, 
    size_t stride, size_t n_surfs, const double *trg, size_t trg_idx_head, 
    size_t trg_idx_tail, double *pot) const
{                  
    PROFILESTART();         
    if(qw != NULL)
        DirectStokesKernel(stride, n_surfs, trg_idx_head, 
            trg_idx_tail, qw, trg, src, den, pot);
    else
        DirectStokesKernel_Noqw(stride, n_surfs, trg_idx_head, 
            trg_idx_tail, trg, src, den, pot);
    
    PROFILEEND("CPUd",0);
    return;
} 

template<enum DeviceType DT>
void Device<DT>::DirectStokes(const float *src, const float *den, const float *qw, 
    size_t stride, size_t n_surfs, const float *trg, size_t trg_idx_head, 
    size_t trg_idx_tail, float *pot) const
{
    PROFILESTART();

#ifdef __SSE2__ 
    if(qw != NULL)
        DirectStokesSSE(stride, n_surfs, trg_idx_head, trg_idx_tail, 
            qw, trg, src, den, pot);
    else
        DirectStokesSSE_Noqw(stride, n_surfs, trg_idx_head, 
            trg_idx_tail, trg, src, den, pot);
    PROFILEEND("SSECPUf",0);
#else
    cerr<<"The SSE instructions in the code are not compatible with this architecture."<<endl;
    PROFILEEND("CPUf",0);
#endif   
}

template<enum DeviceType DT>
template<typename T>
T Device<DT>::Max(T *x_in, size_t length) const
{ 
    PROFILESTART();
    T *max_arr;
    T n_threads;
    
#pragma omp parallel
    {
        if(omp_get_thread_num() == 0)
            max_arr= (T*) this->Malloc(omp_get_num_threads() * sizeof(T));
        
        T max_loc = *x_in;
#pragma omp for
        for(size_t idx = 0;idx<length;idx++)
            max_loc = (max_loc > x_in[idx]) ? max_loc : x_in[idx];
        max_arr[omp_get_thread_num()] =  max_loc;
        
        if(omp_get_thread_num() == 0)
            n_threads = omp_get_num_threads();
    }
    
    T max=max_arr[0];
    for(size_t idx = 0;idx<n_threads;idx++)
        max = (max > max_arr[idx]) ? max : max_arr[idx];
    
    Free(max_arr);

    PROFILEEND("CPU",0);
    return(max);
}

template<enum DeviceType DT>
template<typename T>
T* Device<DT>::Transpose(const T *in, size_t height, size_t width, T *out) const
{ 
    assert(out != in);
    PROFILESTART();
    
#pragma omp parallel for
    for(size_t jj=0;jj<width;++jj)
        for(size_t ii=0;ii<height;++ii)
            out[jj*height + ii] = in[ii * width + jj];

    PROFILEEND("CPU",0);
    return(out);
}

template<enum DeviceType DTlhs, enum DeviceType DTrhs>                         
inline bool operator==(const Device<DTlhs> &lhs, const Device<DTrhs> &rhs)
{
    return((void*) &lhs == (void*) &rhs);
}
