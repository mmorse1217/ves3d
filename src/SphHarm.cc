/**
 * @file   SphHarm.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Jan 31 13:50:30 2010
 * 
 * @brief  The implementation of the SphHarm class. 
 */

#include<iostream>

template<typename scalarType> 
SphHarm<scalarType>::SphHarm()
{
    std::cout<<"Initialization will be here"<<std::endl;
}

template<typename scalarType> 
SphHarm<scalarType>::~SphHarm()
{
    std::cout<<"Destruction will be here"<<std::endl;
}

template <typename scalarType> 
void SphHarm<scalarType>::Derivatives(SHScalars<scalarType> *f_in, 
    SHScalars<scalarType> *Duf_out, SHScalars<scalarType> *Dvf_out, 
    SHScalars<scalarType> *Duuf_out, SHScalars<scalarType> *Duvf_out, 
    SHScalars<scalarType> *Dvvf_out)
{
    std::cout<<"Hello again"<<std::endl;
}

