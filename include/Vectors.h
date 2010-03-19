/**
 * @file   SHVectors.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Thu Jan 28 14:23:48 2010
 * 
 * @brief The header file for the SHVectors class. See
 * src/Vectors.cc for its implementation.
 */

#ifndef _VECTORS_H_
#define _VECTORS_H_

#include "Scalars.h"

//Forward declaration of Vectors 
template <typename T> class Vectors; 

//Function declaration
template<typename T> 
void DotProduct(const Vectors<T>& x_in, 
    const Vectors<T>& y_in, Scalars<T>& xDb_out);

template<typename T> 
void CrossProduct(const Vectors<T>& x_in, 
    const Vectors<T>& y_in, Vectors<T>& xCy_out);

template<typename T> 
void xvpw(const Scalars<T>& x_in, const Vectors<T>& v_in, 
    const Vectors<T>& w_in, Vectors<T>& xvpw_out);

template<typename T> 
void xvpb(const Scalars<T>& a_in, const Vectors<T>& x_in, 
    T b_in, Vectors<T>& xvpb_out); 

template<typename T> 
void uyInv(const Vectors<T>& u_in, const Scalars<T>& y_in, 
    Vectors<T> &uyInv_out); 


/**
 * Three dimensional vectors fields with each component defined on the
 * sphere. The <i>"spherical calculus"</i> is the same as Scalars, but
 * with the added functionality of <i>"vector calculus"</i>. Let
 * <tt>[X_i,Y_i,Z_i]</tt> be a vector field on the sphere, i.e. each
 * component is a scalar field on the sphere. Then the <tt>data_</tt>
 * is organized as <tt>[X_1 Y_1 Z_1 X_2 Y_2 Z_2 ...]</tt>.
 * 
 */
template <typename T> 
class Vectors : public Scalars<T>
{
 public:
    /// Default constructor, no memory allocation.
    Vectors(Device<T> &device_in);
    
    /** 
     * Constructor with size argument. In this case
     * number_of_functions_ will be  3*n_vecs_
     * 
     * @param p_in The grid size variable
     * @param num_vecs_in Number of vectors
     */
    Vectors(Device<T> &device_in, int p_in, int num_vecs_in);
    
    /** 
     * Constructor with size argument and input initializing data.
     * 
     * @param p_in The grid size variable
     * @param num_vecs_in The number of vectors
     * @param vec_data_in The initializing data (assumed to be of correct size).
     */
    Vectors(Device<T> &device_in, int p_in, int num_vecs_in, const T *vec_data_in);

    /** 
     * Returns the length of the vector, that is 3*GetFunLength()
     * 
     * @return The length of each vector 
     */
    int GetVecLength() const;

    void Resize(int n_vecs_in);

    /** 
     * The geometrical dot product on Vectors.
     */
    friend void DotProduct<T>(const Vectors<T> &x_in, 
        const Vectors<T> &y_in, Scalars<T> &xDy_out);
    
    /**
     * The geometrical cross product on Vectors.
     * 
     */
    friend void CrossProduct<T>(const Vectors<T>& x_in, 
        const Vectors<T>& y_in, Vectors<T>& xCy_out);
    
    /**
     * The addition operator for Vectors.
     */
    friend void xvpw<T>(const Scalars<T>& x_in, 
        const Vectors<T>& v_in, const Vectors<T>& w_in, 
        Vectors<T>& xvpw_out); 

    /**
     * The addition operator for Vectors.
     */
    friend void xvpb<T>(const Scalars<T>& x_in, 
        const Vectors<T>& v_in, T b_in, 
        Vectors<T>& xvpb_out); 
    
    friend void uyInv<T>(const Vectors<T>& u_in, 
        const Scalars<T>& y_in, Vectors<T> &uyInv_out); 

  private:
    /// Number of vectors in the class.
    int n_vecs_;
    
    /** 
     * The declaration of the copy constructor, this is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    Vectors(const Vectors<T>& vec_in);
    
    /** 
     * The declaration of the assignment operator, it is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    Vectors<T>& operator=(const Vectors<T>& vec_in);
};

#include "Vectors.cc"
#endif
