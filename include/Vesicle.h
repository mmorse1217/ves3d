/**
 * @file   Vesicle.h
 * @author Rahimian, Abtin <abtin@romario>
 * @date   Tue Feb  2 14:37:48 2010
 * 
 * @brief  The declaration for the class Vesicle.
 */

#ifndef _VESICLE_H_
#define _VESICLE_H_

#include "Surface.h"

template <typename ScalarType> class Vesicle : public Surface<ScalarType>
{
  public:
    ScalarType *kappa_;
    SHScalars<ScalarType> tension_;
    
    Vesicle();
    Vesicle(int p_in, int number_of_vesicles_in);
    Vesicle(int p_in, int number_of_vesicles_in, SHVectors<ScalarType> *x_in);
    ~Vesicle(); 

    void GetBendingForce(SHVectors<ScalarType> *bending_force_out);
    void GetTensileForce(SHVectors<ScalarType> *tensile_force_out);

    //int LiearizedCurvatureOpt(SHVectors<ScalarType> *vec_in, SHScalars<ScalarType> *ininearized_curvature_out); 
    //int LinearizedBendingOpt (SHScalars<ScalarType> *curvature, SHVectors<ScalarType> *bending_force_out);
    //int LinearizedTensionOpt (SHScalars<ScalarType> *tension  , SHVectors<ScalarType> *tensile_force_out);
    
    //friend int StokesMatVec(Vesicle *ves_in,...); //arn?
};

#include "Vesicle.cc"
#endif
