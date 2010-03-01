/**
 * @file   DeviceCpu.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Feb 23 15:28:14 2010
 * 
 * @brief  The implementation of the DeviceCPU class.
 */

//#include "Device.h"

using namespace std;

template<typename T>
T* DeviceCPU<T>::Malloc(size_t size)
{
#ifndef NDEBUG
    cout<<"DeviceCPU::Malloc, size="<<size<<endl;
#endif
    
    return((T*) ::malloc(size));
}

// template<typename T>
// void DeviceCPU<T>::Free(T* ptr)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::Free"<<endl;
// #endif
    
//     ::free(ptr);
// }

// template<typename T>
// T* DeviceCPU<T>::Calloc(size_t num, size_t size)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::Calloc, num="<<num<<", size="<<size<<endl;
// #endif
    
//     return((T*) ::calloc(num,size));
// }

// template<typename T>
// T* DeviceCPU<T>::Memcpy(T* destination, const T* source, size_t num, enum MemcpyKind kind)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::Memcpy"<<endl;
// #endif
 
//     if(destination == source)
//     {
// #ifndef NDEBUG
//         cout<<"  . DeviceCPU::Memcpy, destination == source"<<endl;
// #endif
//         return destination;
//     }
//     else
//     {
//         return((T*) ::memcpy(destination, source, num));
//     }
// }

// template<typename T>  
// T* DeviceCPU<T>::DotProduct(const T* a_in, const T* b_in, int stride, int num_surfs, T* aDb_out)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::DotProduct (single)"<<endl;
// #endif

//     int base, resbase;
//     T dot;

//     for (int surf = 0; surf < num_surfs; surf++) {
//         resbase = surf * stride;
//         base = resbase * 3;
//         for (int s = 0; s < stride; s++) {
//             dot  = a_in[base + s                  ] * b_in[base + s                 ];
//             dot += a_in[base + s + stride         ] * b_in[base + s + stride        ];
//             dot += a_in[base + s + stride + stride] * b_in[base + s +stride + stride];
//             aDb_out[resbase + s] = dot;
//         }
//     }
//     return aDb_out;
// }

// // double* DeviceCPU::DotProduct(const double* a_in, const double* b_in, int stride, int num_surfs, double* aDb_out)
// // {
// // #ifndef NDEBUG
// //     cout<<"DeviceCPU::DotProduct (double)"<<endl;
// // #endif

// //     int base, resbase;
// //     double dot;

// //     for (int surf = 0; surf < num_surfs; surf++) {
// //         resbase = surf * stride;
// //         base = resbase * 3;
// //         for (int s = 0; s < stride; s++) {
// //             dot  = a_in[base + s                  ] * b_in[base + s                 ];
// //             dot += a_in[base + s + stride         ] * b_in[base + s + stride        ];
// //             dot += a_in[base + s + stride + stride] * b_in[base + s +stride + stride];
// //             aDb_out[resbase + s] = dot;
// //         }
// //     }
// //     return aDb_out;
// // }

// template<typename T>
// T* DeviceCPU<T>::CrossProduct(const T* a_in, const T* b_in, int stride, int num_surfs, T* aCb_out)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::CrossProduct (single)"<<endl;
// #endif
    
//     int base, resbase, surf, s;
//     T ax, ay, az, bx, by, bz, x, y, z;

//     for (surf = 0; surf < num_surfs; surf++) {
//         resbase = surf * stride;
//         base = resbase * 3;
//         for (s = 0; s < stride; s++) {
//             ax = a_in[base + s];
//             ay = a_in[base + s + stride];
//             az = a_in[base + s + stride + stride];
//             bx = b_in[base + s];
//             by = b_in[base + s + stride];
//             bz = b_in[base + s + stride + stride];
//             x = ay * bz - az * by;
//             y = az * bx - ax * bz;
//             z = ax * by - ay * bx;
//             aCb_out[base + s] = x;
//             aCb_out[base + s + stride] = y;
//             aCb_out[base + s + stride + stride] = z;
//         }
//     }
//     return aCb_out;
// }

// // double* DeviceCPU::CrossProduct(const double* a_in, const double* b_in, int stride, int num_surfs, double* aCb_out)
// // { 
// // #ifndef NDEBUG
// //     cout<<"DeviceCPU::CrossProduct (double)"<<endl;
// // #endif

// //     int base, resbase, surf, s;
// //     double ax, ay, az, bx, by, bz, x, y, z;

// //     for (surf = 0; surf < num_surfs; surf++) {
// //         resbase = surf * stride;
// //         base = resbase * 3;
// //         for (s = 0; s < stride; s++) {
// //             ax = a_in[base + s];
// //             ay = a_in[base + s + stride];
// //             az = a_in[base + s + stride + stride];
// //             bx = b_in[base + s];
// //             by = b_in[base + s + stride];
// //             bz = b_in[base + s + stride + stride];
// //             x = ay * bz - az * by;
// //             y = az * bx - ax * bz;
// //             z = ax * by - ay * bx;
// //             aCb_out[base + s] = x;
// //             aCb_out[base + s + stride] = y;
// //             aCb_out[base + s + stride + stride] = z;
// //         }
// //     }
// //     return aCb_out;
// // }

// template<typename T>
// T* DeviceCPU<T>::xInv(const T* x_in, int stride, int num_surfs, T* xInv_out)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::xInv (single)"<<endl;
// #endif
//     int length = stride*num_surfs;
//     for (int idx = 0; idx < length; idx++)
//         xInv_out[idx] = 1.0/x_in[idx];
    
//     return xInv_out;
// }

// template<typename T>
// T*  DeviceCPU<T>::AxPy(T a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs , T*  axpy_out)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::AxPy , scalar A (single)"<<endl;
// #endif

//     int length = stride*num_surfs;
//     for (int idx = 0; idx < length; idx++)
//     {
//         axpy_out[idx] = a_in*x_in[idx];
//         axpy_out[idx] += y_in[idx];
//     }
    
//     return axpy_out;
// }

// template<typename T>
// T*  DeviceCPU<T>::AxPy(T a_in, const T*  x_in, T y_in, int stride, int num_surfs , T*  axpy_out)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::AxPy , scalar A, scalar y (single)"<<endl;
// #endif

//     int length = stride*num_surfs;
//     for (int idx = 0; idx < length; idx++)
//     {
//         axpy_out[idx] = a_in*x_in[idx];
//         axpy_out[idx] += y_in;
//     }

//     return axpy_out;
// }

// template<typename T>
// T*  DeviceCPU<T>::AxPy(const T* a_in, const T*  x_in, const T*  y_in, int stride, int num_surfs, T*  axpy_out)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::AxPy , vector A, vector y (single)"<<endl;
// #endif

//     int base, a_base, vec, s, length = 3*stride, idx, a_idx;

//     for (vec = 0; vec < num_surfs; vec++)
//     {
//         base = vec*length;
//         a_base =vec*stride;
        
//         for (s = 0; s < stride; s++) {
            
//             idx = base+s;
//             a_idx = a_base+s;

//             axpy_out[idx]  = a_in[a_idx] * x_in[idx];
//             axpy_out[idx] += y_in[idx];

//             idx +=stride;
//             axpy_out[idx]  = a_in[a_idx] * x_in[idx];
//             axpy_out[idx] += y_in[idx];

//             idx +=stride;
//             axpy_out[idx]  = a_in[a_idx] * x_in[idx];
//             axpy_out[idx] += y_in[idx];
//         }
//     }
//     return axpy_out;
// }

// template<typename T>
// T*  DeviceCPU<T>::AxPy(const T* a_in, const T*  x_in, T y_in, int stride, int num_surfs, T*  axpy_out)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::AxPy , vector A, scalar y(single)"<<endl;
// #endif
//     int base, a_base, vec, s, length = 3*stride, idx, a_idx;

//     for (vec = 0; vec < num_surfs; vec++)
//     {
//         base = vec*length;
//         a_base =vec*stride;
        
//         for (s = 0; s < stride; s++) {
            
//             idx = base+s;
//             a_idx = a_base+s;

//             axpy_out[idx]  = a_in[a_idx] * x_in[idx];
//             axpy_out[idx] += y_in;

//             idx +=stride;
//             axpy_out[idx]  = a_in[a_idx] * x_in[idx];
//             axpy_out[idx] += y_in;

//             idx +=stride;
//             axpy_out[idx]  = a_in[a_idx] * x_in[idx];
//             axpy_out[idx] += y_in;
//         }
//     }

//     return axpy_out;
// }
// template<typename T>
// T* DeviceCPU<T>::xTy(const T* x_in, const T* y_in, int stride, int num_surfs, T* xTy_out)
// {
// #ifndef NDEBUG
//     cout<<"DeviceCPU::xTy (single)"<<endl;
// #endif
    
//     int length = stride*num_surfs, idx;

//     for (idx = 0; idx < length; idx++)
//         xTy_out[idx] = x_in[idx] * y_in[idx];

//     return xTy_out;
// }
