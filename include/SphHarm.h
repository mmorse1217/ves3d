/**
 * @file   SphHarm.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Jan 31 12:51:10 2010
 * 
 * @brief The class declaration for the spherical harmonics transform.
 */

#ifndef _SPHHARM_H_
#define _SPHHARM_H_

#include "SHScalars.h"
#include "cuda_sht.h" 

/** Spherical harmonics transform operator interface.
 * MODIFY THE TEMPLATE TO HAVE THE DEVICE AS A TEMPLATE PARAMETER 
 *
 * COMPOSITE FUNCTIONS TO BE ADDED.
 */
template <typename scalarType> class SphHarm
{  
  private:
    /**
     * 	The grid size variable, equal to the order of spherical
     * 	harmonics expansion of the functions.
     */
    int p_;
    
    /// Number of functions to be operated on.
    int number_of_functions_;
    
    /// The actual transform operator
    cuda_sht  cudaTransClass;
    
    /**
     * The spherical harmonic coefficients of last set of functions
     * passed to the class.  
     */
    scalarType* shc_;
  
 public:
    /// Default constructor.
    SphHarm();
    
    /// Constructor with size argument. THE INTERFACE WITH THE
    /// CUDA/BLAS_SHT IS FLAWED(ADD THE DATA FILE). IT ONLY WORKS FOR P=12.
    SphHarm(int p_in, int number_of_functions_in);

    ///Destructor.
    ~SphHarm();

    /** 
     * Given function set f_in on the sphere, it returns all the first
     * and second derivatives of f_in.
     */
    void AllDerivatives(const SHScalars<scalarType>& f_in, 
        SHScalars<scalarType>& Duf_out, SHScalars<scalarType>& Dvf_out, 
        SHScalars<scalarType>& Duuf_out, SHScalars<scalarType>& Duvf_out, 
        SHScalars<scalarType>& Dvvf_out);

    /** 
     * Given function set f_in on the sphere, it returns all only the
     * derivatives of f_in.
     */
    void FirstDerivatives(const SHScalars<scalarType>& f_in, 
        SHScalars<scalarType>& Duf_out, SHScalars<scalarType>& Dvf_out);

    void Filter(int fiter_freq_in, SHScalars<scalarType>& f_in, SHScalar<scalarType>& filter_f_out);
    void Interpolate(int fiter_freq_in, SHScalars<scalarType>& f_in, SHScalar<scalarType>& filter_f_out);
};

#include "SphHarm.cc"
#endif //_SPHHARM_H_
