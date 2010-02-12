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
#include <cassert>

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
void DotProduct(const SHVectors<ScalarType>& a_in, 
    const SHVectors<ScalarType>& b_in, SHScalars<ScalarType>& aDb_out);

template<typename ScalarType> 
void CrossProduct(const SHVectors<ScalarType>& a_in, 
    const SHVectors<ScalarType>& b_in, SHVectors<ScalarType>& aCb_out);

template<typename ScalarType> 
void AxPy(const SHScalars<ScalarType>& a_in, 
    const SHVectors<ScalarType>& x_in, const SHVectors<ScalarType>& y_in, 
    SHVectors<ScalarType>& c_out); 

template<typename ScalarType> 
void AxPy(const SHScalars<ScalarType>& a_in, 
    const SHVectors<ScalarType>& x_in, ScalarType y_in, 
    SHVectors<ScalarType>& c_out); 

template<typename ScalarType> 
void xDy(const SHVectors<ScalarType>& x_in, const SHScalars<ScalarType>& y_in, 
        SHVectors<ScalarType>& c_out); 

//Class declaration
template <typename ScalarType> 
class SHVectors : public SHScalars<ScalarType>
{
 public:
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
    int GetVecLength() const;

    /** 
     * Sets the private data members p_ and number_of_vectors_ (for
     * now, only before initialization of data_). For the time being
     * it will be impossible to change these values after that the
     * memory for data_ is allocated (it needs an interpolation method
     * to move between to different size).
     */
    void Resize(int p_in, int number_of_vectors_in);

    /** 
     * The vector (geometrical) dot product on SHVectors.
     */
    friend void DotProduct<ScalarType>(const SHVectors<ScalarType>& a_in, 
        const SHVectors<ScalarType>& b_in, SHScalars<ScalarType>& aDb_out);
    
    /**
     * The vector (geometrical) cross product on SHVectors.
     * 
     */
    friend void CrossProduct<ScalarType>(const SHVectors<ScalarType>& a_in, 
        const SHVectors<ScalarType>& b_in, SHVectors<ScalarType>& aCb_out);
    
    /**
     * The addition operator for SHVectors.
     */
    friend void AxPy<ScalarType>(const SHScalars<ScalarType>& a_in, 
        const SHVectors<ScalarType>& x_in, const SHVectors<ScalarType>& y_in, 
        SHVectors<ScalarType>& c_out); 

    /**
     * The addition operator for SHVectors.
     */
    friend void AxPy<ScalarType>(const SHScalars<ScalarType>& a_in, 
        const SHVectors<ScalarType>& x_in, ScalarType y_in, 
        SHVectors<ScalarType>& c_out); 
    
    /**
     * The division (scaling) operator for SHVectors.
     */
    friend void xDy<ScalarType>(const SHVectors<ScalarType>& x_in, 
        const SHScalars<ScalarType>& y_in, SHVectors<ScalarType>& c_out); 

  private:
    /// Number of vectors in the class.
    int number_of_vectors_;
    
    /** 
     * The declaration of the copy constructor, this is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    SHVectors(const SHVectors<ScalarType>& sh_in);
    
    /** 
     * The declaration of the assignment operator, it is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    SHVectors<ScalarType>& operator=(const SHVectors<ScalarType>& sh_in);
};

#include "SHVectors.cc"
#endif
