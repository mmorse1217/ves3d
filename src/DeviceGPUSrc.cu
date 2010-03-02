/**
 * @file   DeviceGPU.cu
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Mon Mar  1 13:09:06 2010
 */

template <typename T>
T* DeviceGPU<T>::Malloc(unsigned long int length)
{
    T* ptr = NULL;
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
    T* ptr = NULL;
    return ptr;
}

template <typename T>
T* DeviceGPU<T>::Memcpy (T* destination, const T* source, unsigned long int num, enum MemcpyKind kind)
{
    return destination;
}

template <typename T>
T* DeviceGPU<T>::DotProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* x_out)
{
    return x_out;
}

template <typename T>
T* DeviceGPU<T>::CrossProduct(const T* u_in, const T* v_in, int stride, int num_surfs, T* w_out)
{
    return w_out;
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
T* DeviceGPU<T>::xvpw(const T* x_in, const T*  v_in, const T*  w_in, int stride, int num_surfs, T*  xvpw_out)
{
    return xvpw_out;
}

template <typename T>
T* DeviceGPU<T>::xvpb(const T* x_in, const T*  v_in, T b_in, int stride, int num_surfs, T*  xvpb_out)
{
    return xvpb_out;
}

// void* DeviceGPU::calloc(size_t num, size_t size)
// {
//     void* ptr = NULL;
//     cudaMalloc(&ptr, size);
//     std::cerr<<"Zero initialization is not implemented."<<std::endl;
//     return(ptr);

// }

// float* DeviceGPU::DotProduct(float* a_in, float* b_in, int stride, int num_vecs, float* aDb_out)
// {
//     return a_in;
// }

// double* DeviceGPU::DotProduct(double* a_in, double* b_in, int stride, int num_vecs, double* aDb_out)
// {
//     std::cerr<<"The double precision operator DotProduct has no implementation for this device."<<std::endl;
//     assert(0);
//     return a_in;
// }

// float* DeviceGPU::CrossProduct(float* a_in, float* b_in, int stride, int num_vecs, float* aDb_out)
// {
//     return a_in;
// }

// double* DeviceGPU::CrossProduct(double* a_in, double* b_in, int stride, int num_vecs, double* aDb_out)
// {
//     std::cerr<<"The double precision operator CrossProduct has no implementation for this device."<<std::endl;
//     assert(0);
//     return a_in;
// }

// float* DeviceGPU::AxPy(float* a_in, float* x_in, float* y_in, int stride, int num_vecs, float* axpy_out)
// {
//     return a_in;
// }

// double* DeviceGPU::AxPy(double* a_in, double* x_in, double* y_in, int stride, int num_vecs, double* axpy_out)
// {
//     std::cerr<<"The double precision operator AxPy has no implementation for this device."<<std::endl;
//     assert(0);
//     return a_in;
// }



// // //cudaMemcpy(a_dp, a, sizeof(real) * vecSize * num_surfs, cudaMemcpyHostToDevice);
// // //cudaMemcpy(aDb, aDb_dp, sizeof(real) * stride * num_surfs, cudaMemcpyDeviceToHost);

// // void* DeviceGpu::calloc(size_t num, size_t size)
// // {
// //     return(::calloc(num,size));
// // }
