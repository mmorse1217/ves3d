/**
 * @file   Scalars.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 26 14:44:00 2010
 * 
 * @brief The header file for the Scalars class. The implementation
 * of this class is in src/Scalars.cc. 
 */

#ifndef _SCALARS_H_
#define _SCALARS_H_

#include "Device.h"
#include <cassert>

//Froward declaration
template<typename T> class Scalars;

//Friend functions' declaration
template<typename T> 
void axpy(T a_in, const Scalars<T> &x_in, const Scalars<T> &y_in, Scalars<T> &axpy_out);

template<typename T> 
void axpb(T a_in, const Scalars<T> &x_in, T y_in, Scalars<T> &axpb_out);

template<typename T> 
void xy(const Scalars<T> &x_in, const Scalars<T> &y_in, Scalars<T> &xy_out);

template<typename T> 
void xyInv(const Scalars<T> &x_in, const Scalars<T> &y_in, Scalars<T>  &xyInv_out);

/**
 * @brief Scalars is a class of scalar fields/functions defined on a
 * Gauss-uniform grid on a sphere.
 *
 * An instance of this class may hold multiple scalar fields so that
 * they can be streamed for any manipulation.
 * 
 * The implementation is in src/Scalars.cc and the tester is
 * test/ScalarsTest.cc
 */
template <typename T> 
class Scalars
{
  public:
    /// The reference to the device that the data member will reside.
    Device<T> &device_;
    
    /**
     * The grid size variable, equal to the order of spherical
     * harmonics expansion of the functions.
     */
    int p_; 
  
    /// Number of functions in the class. 
    int n_funs_;   

    /// The array that holds the function values.
    T *data_;          
    
    /// The maximum number of functions that fit in the allocated
    /// memory. This may as well be different from n_funs_.
    int max_n_funs_;
    
  private:
    /** 
     * The declaration of the copy constructor, this is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    Scalars(const Scalars<T> &sc_in);
    
    /** 
     * The declaration of the assignment operator, it is declared as
     * private so as to disallow any passing by value to
     * functions. There will be no implementation for this method.
     */
    Scalars<T>& operator=(const Scalars<T> &sc_in);


  public:
    //No default constructor, because the device should be known.
    
    ///The device is set but no memory is allocated.
    Scalars(Device<T> &device_in);
    
    /// Constructor with size arguments, the memory is allocated but not
    /// initialized.
    Scalars(Device<T> &device_in, int p_in, int n_funs_in);

    /// Constructor with size arguments and an input data to
    /// initialize the data stored in the class. <b>The input
    /// <code>data_in</code> is assumed to be sitting on the host
    /// (with respect to the device, of course)</b>.
    Scalars(Device<T> &device_in, int p_in, int n_funs_in, const T *data_in); 

    /// Default destructor, the memory is freed.
    virtual ~Scalars();
    
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
     * GetDataLength(). . <b>The input <code>data_in</code> is assumed
     * to be sitting on the host</b>.
     */
    void SetData(const T *data_in); 
    
    /** 
     * Set the values of a function at a given index. fun_in is
     * assumed to be of length that of returned by
     * GetFunLength(). <b>The input <code>fun_in</code> is assumed to
     * be sitting on the host</b>.
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

    T Max();

    /** 
     * Sets the private data members p_ and n_funs_ (for now, only
     * before initialization of data_). For the time being it will be
     * impossible to change these values after that the memory for
     * data_ is allocated (it needs an interpolation method to move
     * between to different size). It also allocates memory for data_.
     */
    virtual void Resize(int n_funs_in);

    friend void axpy<T>(T a_in, const Scalars<T> &x_in, 
        const Scalars<T> &y_in, Scalars<T> &axpy_out);

    friend void axpb<T>(T a_in, const Scalars<T> &x_in, T b_in, 
        Scalars<T> &axpb_out);

    friend void xy<T>(const Scalars<T> &x_in, const Scalars<T> &y_in, 
        Scalars<T> &xy_out);

    friend void xyInv<T>(const Scalars<T> &x_in, const Scalars<T> &y_in,
        Scalars<T> &xyInv_out);

};

#include "Scalars.cc"

#endif //_SCALARS_H_
