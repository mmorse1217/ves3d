/**
 * @file   SHVectors.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Fri Jan 29 09:36:58 2010
 * 
 * @brief  The implementation of SHVectors.
 */

#include <stdexcept>
#include<iostream>

// Constructors
template<typename ScalarType> 
SHVectors<ScalarType>::SHVectors() :
    number_of_vectors_(0){}

template<typename ScalarType> 
SHVectors<ScalarType>::SHVectors(int p_in, int num_vecs_in) :
    SHScalars<ScalarType>( p_in, 3*num_vecs_in) , number_of_vectors_(num_vecs_in){}

template<typename ScalarType> 
SHVectors<ScalarType>::SHVectors(int p_in, int num_vecs_in, const ScalarType *vec_data_in) :
    SHScalars<ScalarType>( p_in, 3*num_vecs_in, vec_data_in) , number_of_vectors_(num_vecs_in){}

// Utility functions
template<typename ScalarType> 
int SHVectors<ScalarType>::GetVecLength()
{
    return(3*SHScalars<ScalarType>::GetFunLength());
}


template<typename ScalarType> 
void DotProduct(SHVectors<ScalarType> *a_in, 
    SHVectors<ScalarType> *b_in, SHScalars<ScalarType> *aDb_out)
{

    std::cout<<" The dot product is not optimized!!"<<std::endl;

    if(a_in->GetDataLength() != b_in->GetDataLength())
        throw std::length_error(" Input vector sizes do not match.");

    if(a_in->GetDataLength() > 3*aDb_out->GetDataLength())
        throw std::length_error(" Output vector is too small.");
    
    int funLen = a_in->GetFunLength();
    int numVec = a_in->number_of_vectors_;
    int idx;
    const ScalarType *a_data = a_in->GetData();
    const ScalarType *b_data = a_in->GetData();
    ScalarType *dot = new ScalarType[numVec*funLen];
    
    for(int ii=0;ii<numVec;++ii)
    {
        a_data = a_in->GetFunctionAt(3*ii);
        b_data = b_in->GetFunctionAt(3*ii);

        for(int jj=0;jj<funLen;++jj)
        {
            idx = ii*funLen +jj;
            
            dot[idx] = a_data[jj         ]*b_data[jj         ];
            dot[idx]+= a_data[jj+  funLen]*b_data[jj+  funLen];
            dot[idx]+= a_data[jj+2*funLen]*b_data[jj+2*funLen];
            
        }
    }
    aDb_out->SetData(dot);

}
template<typename ScalarType> 
void CrossProduct(SHVectors<ScalarType> *a_in, 
    SHVectors<ScalarType> *b_in, SHVectors<ScalarType> *aCb_out)
{
    std::cout<<" The cross product is not optimized!!"<<std::endl;

    if(a_in->GetDataLength() != b_in->GetDataLength() ||
        b_in->GetDataLength() !=aCb_out->GetDataLength())
        throw std::length_error(" Input vector sizes do not match.");
    
    int funLen = a_in->GetFunLength();
    int numVec = a_in->number_of_vectors_;
    int idx;
    const ScalarType *a_data = a_in->GetData();
    const ScalarType *b_data = a_in->GetData();
    ScalarType *cross = new ScalarType[aCb_out->GetDataLength()];
    
    for(int ii=0;ii<numVec;++ii)
    {
        a_data = a_in->GetFunctionAt(3*ii);
        b_data = b_in->GetFunctionAt(3*ii);

        for(int jj=0;jj<funLen;++jj)
        {
            idx = ii*funLen +jj;
            
            cross[idx         ] = a_data[jj+  funLen]*b_data[jj+2*funLen] - a_data[jj+2*funLen]*b_data[jj+  funLen];
            cross[idx+  funLen] = a_data[jj+2*funLen]*b_data[jj         ] - a_data[jj         ]*b_data[jj+2*funLen];
            cross[idx+2*funLen] = a_data[jj         ]*b_data[jj+  funLen] - a_data[jj+  funLen]*b_data[jj         ];
            
        }
    }
    aCb_out->SetData(cross);
    
}

