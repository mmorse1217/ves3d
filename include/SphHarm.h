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
/// Spherical harmonics
template <typename scalarType> class SphHarm
{  
 private:
 public:
    SphHarm();
    //SphHarm(int p_in, int number_of_functions_in);
    ~SphHarm();

    void Derivatives(SHScalars<scalarType> *f_in, 
        SHScalars<scalarType> *Duf_out, SHScalars<scalarType> *Dvf_out, 
        SHScalars<scalarType> *Duuf_out, SHScalars<scalarType> *Duvf_out, 
        SHScalars<scalarType> *Dvvf_out);
    //int PointsToFrequency(SHScalars *f_in, SHScalars *f_out);
    //int FrequencyToPoints(SHScalars *f_in, SHScalars *f_out);
    //int Resample(SHScalars *f_in, SHScalars *f_out);
    //int ResampleWithScaling(SHScalars *f_in, int *p_scaling_in, SHScalars *f_out);
    
    // composite functions
};

#include "SphHarm.cc"
#endif
