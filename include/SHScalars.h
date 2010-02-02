/**
 * @file   SHScalars.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 26 14:44:00 2010
 * 
 * @brief The headier file for the SHScalars class. The implementation
 * of this class is in src/SHScalars.cc. 
 */

#ifndef _SHSCALARS_H_
#define _SHSCALARS_H_

#include <stdexcept>

/**
 * @brief SHScalars is a class of scalar field/functions defined on a
 * Gauss-uniform grid on a sphere.
 *
 * An instance of this class may hold multiple scalar fields so that
 * they can be streamed for any manipulation, such as differentiation.
 * 
 * The implementation is in src/SHScalars.cc and the tester is
 * test/SHScalarsTest.cc
 */
template <typename ScalarType> 
class SHScalars
{
  protected:
  public:
    /// The grid size variable, equal to the order of spherical harmonics expansion of the functions.
    int p_; 
  
    /// Number of functions in the class. 
    int number_of_functions_;   

    /// The array that holds the function values.
    ScalarType *data_;          
  
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
    SHScalars(int p_in, int num_funs_in, const ScalarType *data_in); 

    /// Default deconstructor, the memory is freed.
    virtual ~SHScalars();
    
    /** 
     * Calculates the length of each field as a function of p_. 
     *
     * @return The length of each field.
     */
    int GetFunLength();
    
    /** 
     * Calculates the length of data_ array.
     *
     * @return The length of data_ array.
     */
    int GetDataLength();

    /** 
     * The setter function for data_.
     * 
     * @param data_in pointer to the head of input array. It is
     * assumed to be of length at least that of returned by
     * GetDataLength().
     */
    void SetData(const ScalarType *data_in); 
    
    /** 
     * Set the values of a function at a given index. fun_in is
     * assumed to be of length that of returned by GetFunLength().
     * 
     * @param fun_in The array of function values.
     * @param fun_idx The index of the target location, starting from 0.
     * 
     */
    void SetFunctionAt(const ScalarType *fun_in, int fun_idx_in);
    
    /** 
     * The getter function for a field at a given index.
     * 
     * @param fun_idx_in The index of the function, starts from 0.
     * @return A constant pointer to the Function.
     */
    const ScalarType* GetFunctionAt(int fun_idx_in);
};

#include "SHScalars.cc"
#endif
