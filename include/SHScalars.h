/**
 * @file   SHScalars.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 26 14:44:00 2010
 * 
 * @brief The headier file for the SHScalars class. The implementation
 * of this class is in src/SHScalars.cc. 
 */

#ifndef _SCALARS_H_
#define _SCALARS_H_

#include <cassert>
#include "Device.h"

//Forward declaration
template <typename T> class SHScalars; 

//Function declaration
template<typename T> 
void AxPy(T a_in, const SHScalars<T>& x_in, 
    const SHScalars<T>& y_in, SHScalars<T>& res_out); 

template<typename T> 
void AxPy(T a_in, const SHScalars<T>& x_in, 
    T y_in, SHScalars<T>& res_out);

template<typename T> 
void xTy(const SHScalars<T>& x_in, const SHScalars<T>& y_in, 
    SHScalars<T>& xTy_out);

template<typename T> 
void xDy(const SHScalars<T>& x_in, const SHScalars<T>& y_in, 
    SHScalars<T>& xDy_out);

/**
 * @brief SHScalars is a class of scalar fields/functions defined on a
 * Gauss-uniform grid on a sphere.
 *
 * An instance of this class may hold multiple scalar fields so that
 * they can be streamed for any manipulation, such as differentiation.
 * 
 * The implementation is in src/SHScalars.cc and the tester is
 * test/SHScalarsTest.cc
 */
template <typename T> 
class SHScalars
{
  protected:
  public:
    /// The array that holds the function values.
    T *data_;          
  
    /// Default constructor, no memory is allocated.
    SHScalars();                

    /** 
     * Constructor with size arguments, the memory is allocated but not
     * initialized.
     * 
     * @param p_in The grid size variable
     * @param num_funs_in Number of functions
     */
    SHScalars(int p_in, int num_funs_in);

    /** 
     * Constructor with size arguments and an input data to initialize
     * the data stored in the class.
     * 
     * @param p_in The grid size variable
     * @param num_funs_in Number of functions
     * @param data_in The data array to initialize the data in the
     * file. data_in is assumed to be of length at least that of
     * returned by GetDataLength().
     */
    SHScalars(int p_in, int num_funs_in, const T *data_in); 

    /// Default destructor, the memory is freed.
    virtual ~SHScalars();
    
    /** 
     * Calculates the length of each field as a function of p_. 
     *
     * @return The length of each field.
     */
    int GetFunLength() const;
    
    /** 
     * Calculates the length of data_ array.
     *
     * @return The length of data_ array.
     */
    int GetDataLength() const;

    /** 
     * The setter function for data_.
     * 
     * @param data_in pointer to the head of input array. It is
     * assumed to be of length at least that of returned by
     * GetDataLength().
     */
    void SetData(const T *data_in); 
    
    /** 
     * Set the values of a function at a given index. fun_in is
     * assumed to be of length that of returned by GetFunLength().
     * 
     * @param fun_in The array of function values.
     * @param fun_idx The index of the target location, starting from 0.
     * 
     */
    void SetFunctionAt(const T *fun_in, int fun_idx_in);
    
    /** 
     * The getter function for a field at a given index.
     * 
     * @param fun_idx_in The index of the function, starts from 0.
     * @return A constant pointer to the Function.
     */
    const T* GetFunctionAt(int fun_idx_in) const;


    /** 
     * The square root operator for the class. It calculates and saves
     * the point wise square roots.
     */
    void Sqrt();

    /** 
     * Sets the private data members p_ and number_of_functions_ (for
     * now, only before initialization of data_). For the time being
     * it will be impossible to change these values after that the
     * memory for data_ is allocated (it needs an interpolation method
     * to move between to different size). It also allocates memory
     * for data_.
     */
    virtual void Resize(int p_in, int number_of_functions_in);

    /** 
     * The (overloaded) addition operator for the class type
     * SHScalars. There will be no size check inside of the function,
     * and the input is assumed to be of the correct size. The
     * function calculates c = a*x + y.
     */
    friend void AxPy<T>(T a_in, 
        const SHScalars<T>& x_in, const SHScalars<T>& y_in, 
        SHScalars<T>& res_out);

    /** 
     * The (overloaded) addition operator for the class type
     * SHScalars. There will be no size check inside of the function,
     * and the input is assumed to be of the correct size. The
     * function calculates c = a*x + y.
     */
    friend void AxPy<T>(T a_in, 
        const SHScalars<T>& x_in, T y_in, 
        SHScalars<T>& res_out);

    /** 
     * The multiplication operator for the class type SHScalars. There
     * will be no size check inside of the function, and the input is
     * assumed to be of the correct size. The function calculates c =
     * x*y (pointwise).
     */
    friend void xTy<T>(const SHScalars<T>& x_in, 
        const SHScalars<T>& y_in, SHScalars<T>& xTy_out);

    /** 
     * The division operator for the class type SHScalars. There will
     * be no size check inside of the function, and the input is
     * assumed to be of the correct size. The function calculates c =
     * x/y (pointwise).
     */
    friend void xDy<T>(const SHScalars<T>& x_in, 
        const SHScalars<T>& y_in, SHScalars<T>& xDy_out);


  private:
    /**
     * The grid size variable, equal to the order of spherical
     * harmonics expansion of the functions.
     */
    int p_; 
  
    /// Number of functions in the class. 
    int number_of_functions_;   

    /** 
     * The declaration of the copy constructor, this is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    SHScalars(const SHScalars<T>& sh_in);
    
    /** 
     * The declaration of the assignment operator, it is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    SHScalars<T>& operator=(const SHScalars<T>& sh_in);

    /**
     * Allocated the memory in the heap for the data_ member. The size
     * of the allocation will be equal to the value returned by
     * GetDataLength().
     */
    void AllocateMemory();
};

#include "SHScalars.cc"

#endif //_SHSCALARS_H_
