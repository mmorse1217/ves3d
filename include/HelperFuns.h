#ifndef _HELPERFUNS_H_
#define _HELPERFUNS_H_

#include "Device.h"
#include "Scalars.h"
#include "Vectors.h"
#include <iostream>

template<typename T, enum DeviceType DT>
inline bool AreCompatible(const Scalars<T,DT> &lhs, const Scalars<T,DT> &rhs);

template<typename T, enum DeviceType DT>
inline bool AreCompatible(const Vectors<T,DT> &lhs, const Scalars<T,DT> &rhs);

template<typename T, enum DeviceType DT>
inline bool AreCompatible(const Scalars<T,DT> &lhs, const Vectors<T,DT> &rhs);

template<typename T, enum DeviceType DT>
inline bool AreCompatible(const Vectors<T,DT> &lhs, const Vectors<T,DT> &rhs);

template<typename T, enum DeviceType DT>
inline void Sqrt(const Scalars<T,DT> &x_in, Scalars<T,DT> &sqrt_out);

template<typename T, enum DeviceType DT>
inline void xy(const Scalars<T,DT> &x_in, 
               const Scalars<T,DT> &y_in, 
                     Scalars<T,DT> &xy_out);

template<typename T, enum DeviceType DT>
inline void xyInv(const Scalars<T,DT> &x_in, 
                  const Scalars<T,DT> &y_in, 
                        Scalars<T,DT> &xyInv_out);

template<typename T, enum DeviceType DT>
inline void xInv(const Scalars<T,DT> &x_in, 
                       Scalars<T,DT> &xInv_out);
    
template<typename T, enum DeviceType DT>
inline void axpy(T a_in, const Scalars<T,DT> &x_in, 
                         const Scalars<T,DT> &y_in, 
                               Scalars<T,DT> &axpy_out);

template<typename T, enum DeviceType DT>
inline void axpy(T a_in, const Scalars<T,DT> &x_in, 
                               Scalars<T,DT> &axpy_out);

template<typename T, enum DeviceType DT>
inline void Reduce(const Scalars<T,DT> &x_in, 
                   const Scalars<T,DT> &w_in, 
                   const Scalars<T,DT> &quad_w_in, 
                         Scalars<T,DT> &int_x_dw);
    
template<typename T, enum DeviceType DT>
inline void Reduce(const Scalars<T,DT> &w_in, 
                   const Scalars<T,DT> &quad_w_in, 
                         Scalars<T,DT> &int_x_dw);

template<typename T, enum DeviceType DT>
inline T Max(const Scalars<T,DT> &x_in);

template<typename T, enum DeviceType DT>
inline void DotProduct(const Vectors<T,DT> &u_in, 
                       const Vectors<T,DT> &v_in, 
                             Scalars<T,DT> &x_out);

template<typename T, enum DeviceType DT>
inline void CrossProduct(const Vectors<T,DT> &u_in, 
                         const Vectors<T,DT> &v_in, 
                               Vectors<T,DT> &w_out);

template<typename T, enum DeviceType DT>
inline void uyInv(const Vectors<T,DT> &u_in, 
                  const Scalars<T,DT> &y_in, 
                        Vectors<T,DT> &uyInv_out);
    
template<typename T, enum DeviceType DT>
inline void avpw(const Scalars<T,DT> &a_in, 
                 const Vectors<T,DT> &v_in, 
                 const Vectors<T,DT> &w_in, 
                       Vectors<T,DT> &avpw_out);
 
template<typename T, enum DeviceType DT>
inline void xvpw(const Scalars<T,DT> &x_in, 
                 const Vectors<T,DT> &v_in,
                 const Vectors<T,DT> &w_in, 
                       Vectors<T,DT> &xvpw_out);

template<typename T, enum DeviceType DT>
inline void xv(const Scalars<T,DT> &x_in, 
               const Vectors<T,DT> &v_in,
                     Vectors<T,DT> &xvpw_out);
    
template<typename T, enum DeviceType DT>
inline void Reduce(const Vectors<T,DT> &x_in, 
                   const Scalars<T,DT> &w_in, 
                   const Scalars<T,DT> &quad_w_in, 
                         Vectors<T,DT> &int_x_dw);
    
template<typename T, enum DeviceType DT>
inline void ShufflePoints(const Vectors<T,DT> &x_in, 
                                Vectors<T,DT> &x_out);

template<typename T, enum DeviceType DT>
std::ostream& operator<<(std::ostream& output, Scalars<T,DT> &sc);

#include "HelperFuns.cc"

#endif _HELPERFUNS_H_
