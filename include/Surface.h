/**
 * @file   Surface.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:23:11 2010
 * 
 * @brief  The declaration for the Surface class.
 */

#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "HelperFuns.h"
#include "SHTrans.h"
#include "SurfaceParams.h"
///@todo These two headers need to me moved out
#include "DataIO.h"
#include "OperatorsMats.h"

template <typename ScalarContainer, typename VectorContainer> 
class Surface
{
  public:
    typedef typename ScalarContainer::value_type value_type;
    typedef VectorContainer Vec;
    typedef ScalarContainer Sca;

  private:
    SurfaceParams<value_type> params_;
    SHTrans<value_type, CPU> sht_;
    size_t capacity_;
  
  public:
    Surface(SurfaceParams<value_type> &params_in, 
        const OperatorsMats<value_type> &mats);
    ~Surface();

    VectorContainer x_;
    VectorContainer normal_;

    ScalarContainer w_;
    ScalarContainer h_;
    ScalarContainer k_;
   
    void UpdateFirstForms();
    void UpdateAll();

    void Grad(const ScalarContainer &f_in, VectorContainer &grad_f_out) const;
    void Div(const VectorContainer &f_in, ScalarContainer &div_f_out) const;
 
    void Resize(int n_surfs_in);

    ///@todo The size function should be added.
  private:
    VectorContainer cu_;
    VectorContainer cv_;
 
    ///@todo remove the work vectors
    mutable value_type *shc, *work_arr;    
    mutable ScalarContainer S1, S2, S3, S4, S5, S6;
    mutable VectorContainer V1, V2;
    //ScalarContainer S10;
    //VectorContainer V10, V11, V12, V13;    
    
    //void Reparam();
    //T *alpha_p;
    //T *all_rot_mats_;
    //T *quad_weights_;
    //void StokesMatVec(const VectorContainer &density_in, VectorContainer &velocity_out);
    //void GetTension(const VectorContainer &v_in, const VectorContainer &v_ten_in, T *tension_out);
    //T Area();
    //void Volume();
    //void Populate(const T *centers);
    //bool IsAccurate();
    //T* GetCenters(T* cnts);
    //void UpdateNormal();
    //T max_init_area_;  
    //T *alpha_q;
    //Rotation
    //ScalarContainer w_sph_;
    //T *rot_mat;
    //T *sing_quad_weights_;
    //double StokesMatVec_time_;
};

#include "Surface.cc"

#endif //_SURFACE_H_
