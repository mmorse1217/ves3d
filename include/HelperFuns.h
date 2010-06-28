#ifndef _HELPERFUNS_H_
#define _HELPERFUNS_H_

#include <cassert>
#include <iostream>
#include "Device.h"
#include "SHTrans.h"

namespace containers 
{
    template<typename lhsContainer, typename rhsContainer>
    inline bool AreCompatible(const lhsContainer &lhs,
        const rhsContainer &rhs);

    template<typename ScalarContainer>
    inline void Sqrt(const ScalarContainer &x_in, 
        ScalarContainer &sqrt_out);

    template<typename ScalarContainer>
    inline void xy(const ScalarContainer &x_in, 
        const ScalarContainer &y_in, ScalarContainer &xy_out);

    template<typename ScalarContainer>
    inline void xyInv(const ScalarContainer &x_in, 
        const ScalarContainer &y_in, ScalarContainer &xyInv_out);

    template<typename ScalarContainer>
    inline void xInv(const ScalarContainer &x_in, 
        ScalarContainer &xInv_out);
    
    template<typename ScalarContainer>
    inline void axpy(typename ScalarContainer::value_type a_in, 
        const ScalarContainer &x_in, const ScalarContainer &y_in, 
        ScalarContainer &axpy_out);

    template<typename ScalarContainer>
    inline void axpy(typename ScalarContainer::value_type a_in, 
        const ScalarContainer &x_in, ScalarContainer &axpy_out);

    template<typename ScalarContainer>
    inline typename ScalarContainer::value_type Max(const ScalarContainer &x_in);

    template<typename ScalarContainer, typename VectorContainer>
    inline void GeometricDot(const VectorContainer &u_in, 
        const VectorContainer &v_in, ScalarContainer &x_out);

    template<typename VectorContainer>
    inline void GeometricCross(const VectorContainer &u_in, 
        const VectorContainer &v_in, VectorContainer &w_out);

    template<typename ScalarContainer, typename VectorContainer>
    inline void uyInv(const VectorContainer &u_in, 
        const ScalarContainer &y_in, VectorContainer &uyInv_out);
    
    template<typename ScalarContainer, typename VectorContainer>
    inline void avpw(const ScalarContainer &a_in, 
        const VectorContainer &v_in, const VectorContainer &w_in, 
        VectorContainer &avpw_out);
 
    template<typename ScalarContainer, typename VectorContainer>
    inline void xvpw(const ScalarContainer &x_in, 
        const VectorContainer &v_in, const VectorContainer &w_in, 
        VectorContainer &xvpw_out);

    template<typename ScalarContainer, typename VectorContainer>
    inline void xv(const ScalarContainer &x_in, 
        const VectorContainer &v_in, VectorContainer &xvpw_out);
    
    template<typename ScalarContainer, typename IntegrandContainer>
    inline void Reduce(const IntegrandContainer &x_in, 
        const ScalarContainer &w_in, const ScalarContainer &quad_w_in, 
        IntegrandContainer &x_dw);

    template<typename Container>
    inline void Reduce(const Container &w_in, const Container &quad_w_in, 
        Container &dw);
    
    template<typename VectorContainer>
    inline void ShufflePoints(const VectorContainer &x_in, 
        VectorContainer &x_out);

    template<typename ScalarContainer>
    inline void CircShift(const typename ScalarContainer::value_type *x_in,
        int shift, ScalarContainer &x_out);
    
    template<typename VectorContainer>
    inline void Populate(VectorContainer &x,VectorContainer &centers);

    template<typename Container>
    inline typename Container::value_type MaxAbs(Container &x);

    ///@todo this need to be implemented for the GPU
    template<typename ScalarContainer>
    typename ScalarContainer::value_type AlgebraicDot( 
        const ScalarContainer &x, const ScalarContainer &y);

    template<typename Container, typename SHT>
    inline void Resample(const Container &xp, const SHT &shtp, const SHT &shtq, 
        Container &shcpq, Container &wrkpq, Container &xq);
    
#include "HelperFuns.cc"
}

#endif _HELPERFUNS_H_
