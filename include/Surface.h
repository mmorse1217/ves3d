/**
 * @file   Surface.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:23:11 2010
 * 
 * @brief  The declaration for the Surface class.
 */

#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "SHScalars.h"
#include "SHVectors.h"
#include "SphHarm.h"

template <typename ScalarType> class Surface
{
  private:
    ScalarType *area_, *volume_;
  public:
    int p_;
    int number_of_surfs_;
    SHVectors<ScalarType> x_, normal_;
    SHScalars<ScalarType> h_, w_, k_, cu_, cv_;

    Surface();
    Surface(int p_in, int number_of_surfs_in);
    Surface(int p_in, int number_of_surfs_in, const SHVectors<ScalarType> &x_in);

    void SetX(const SHVectors<ScalarType> &x_in);

    //~Surface();
    //int SurfGrad(SHScalars *scalar_in, SHVectors* surf_grad_vec_out);
    //int SurdDiv(SHVectors* vector_in, SHScalars* surf_div_out);
};

#include "Surface.cc"
#endif
