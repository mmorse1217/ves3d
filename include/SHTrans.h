/**
 * @file   SHTrans.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Jan 31 12:51:10 2010
 * 
 * @brief The class declaration for the spherical harmonics transform.
 */

#ifndef _SHTRANS_H_
#define _SHTRANS_H_

#include "Device.h"
#include "Scalars.h"
#include "cuda_sht.h" 

/** Spherical harmonics transform operator interface.
 *
 * @todo composite functions to be added.
 */
template <typename T> class SHTrans
{  
  public:
    /**
     * 	The grid size variable, equal to the order of spherical
     * 	harmonics expansion of the functions.
     */
    int p_;
    
    /// Number of functions to be operated on.
    int n_funs_;
    
    /// The actual transform operator
    cuda_sht  cudaTransClass;
    
    /**
     * The spherical harmonic coefficients of last set of functions
     * passed to the class.  
     */
    T* shc_;
  
    Device<T> &device_;

 public:
    /// Default constructor.
    SHTrans(Device<T> &device_in);
    
    /// Constructor with size argument. THE INTERFACE WITH THE
    /// @todo CUDA/BLAS_SHT is flawed(add the data file). it only works for p=12.
    SHTrans(Device<T> &device_in, int p_in, int n_funs_in);

    ///Destructor.
    ~SHTrans();

    /** 
     * Given function set f_in on the sphere, it returns all the first
     * and second derivatives of f_in.
     */
    void AllDerivatives(const Scalars<T> &f_in, 
        Scalars<T> &Duf_out, Scalars<T> &Dvf_out, 
        Scalars<T> &Duuf_out, Scalars<T> &Duvf_out, 
        Scalars<T> &Dvvf_out);

    /** 
     * Given function f_in on the sphere, it returns only the first
     * derivatives of f_in.
     */
    void FirstDerivatives(const Scalars<T> &f_in, 
        Scalars<T> &Duf_out, Scalars<T> &Dvf_out);

    void Filter(int fiter_freq_in, Scalars<T> &f_in, Scalars<T>  &filtered_f_out);
    void Interpolate(int fiter_freq_in, Scalars<T> &f_in, Scalars<T> &filtered_f_out);
};

#include "SHTrans.cc"
#endif //_SHTRANS_H_
