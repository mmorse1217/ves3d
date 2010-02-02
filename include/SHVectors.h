/**
 * @file   SHVectors.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Thu Jan 28 14:23:48 2010
 * 
 * @brief The header file for the SHVectors class. See
 * src/SHVectors.cc for its implementation.
 */

#ifndef _SHVECTORS_H_
#define _SHVECTORS_H_

#include "SHScalars.h"
#include <stdexcept>
#include<iostream>

/**
 * Three dimensional vectors fields with each component defined on the
 * sphere. The "spherical calculus" is the same as SHScalars, but with
 * the added functionality of "vector calculus". Let [X_i,Y_i,Z_i] be
 * a vector field on the sphere, i.e. each component is a scalar field
 * on the sphere. Then the data_ is organized as [X_1 Y_1 Z_1 X_2 Y_2
 * Z_2 ...].
 * 
 */

//Forward declaration of SHVectors 
template <typename ScalarType> class SHVectors; 

//Function declaration
template<typename ScalarType> 
void DotProduct(SHVectors<ScalarType> *a_in, 
    SHVectors<ScalarType> *b_in, SHScalars<ScalarType> *aDb_out);

template<typename ScalarType> 
void CrossProduct(SHVectors<ScalarType> *a_in, 
    SHVectors<ScalarType> *b_in, SHVectors<ScalarType> *aCb_out);

template<typename ScalarType> 
void AxPy(SHScalars<ScalarType> *a_in, 
        SHVectors<ScalarType> *x_in, SHVectors<ScalarType> *y_in, 
        SHVectors<ScalarType> *c_out); 

template<typename ScalarType> 
void AxPy(ScalarType a_in, 
        SHVectors<ScalarType> *x_in, SHVectors<ScalarType> *y_in, 
        SHVectors<ScalarType> *c_out); 

template<typename ScalarType> 
void AxPy(SHScalars<ScalarType> *a_in, 
        SHVectors<ScalarType> *x_in, ScalarType y_in, 
        SHVectors<ScalarType> *c_out); 

template<typename ScalarType> 
void AxPy(ScalarType a_in, 
        SHVectors<ScalarType> *x_in, ScalarType y_in, 
        SHVectors<ScalarType> *c_out); 


//Class declaration
template <typename ScalarType> 
class SHVectors : public SHScalars<ScalarType>
{
 public:
    /// Number of vectors in the class.
    int number_of_vectors_;

    /// Default constructor, no memory allocation.
    SHVectors();
    
    /** 
     * Constructor with size argument. In this case
     * number_of_functions_ will be  3*number_of_vectors_
     * 
     * @param p_in The grid size variable
     * @param num_vecs_in Number of vectors
     */
    SHVectors(int p_in, int num_vecs_in);
    
    /** 
     * Constructor with size argument and input initializing data.
     * 
     * @param p_in The grid size variable
     * @param num_vecs_in The number of vectors
     * @param vec_data_in The initializing data (assumed to be of correct size).
     */
    SHVectors(int p_in, int num_vecs_in, const ScalarType *vec_data_in);

    /** 
     * Returns the length of the vector, that is 3*GetFunLength()
     * 
     * @return The length of each vector 
     */
    int GetVecLength();

    friend void DotProduct<ScalarType>(SHVectors<ScalarType> *a_in, 
        SHVectors<ScalarType> *b_in, SHScalars<ScalarType> *aDb_out);
    
    friend void CrossProduct<ScalarType>(SHVectors<ScalarType> *a_in, 
        SHVectors<ScalarType> *b_in, SHVectors<ScalarType> *aCb_out);
    
    friend void AxPy<ScalarType>(SHScalars<ScalarType> *a_in, 
        SHVectors<ScalarType> *x_in, SHVectors<ScalarType> *y_in, 
        SHVectors<ScalarType> *c_out); 
};

#include "SHVectors.cc"
#endif
