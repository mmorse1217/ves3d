/**
 * @file   Surface.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:23:11 2010
 * 
 * @brief  The declaration for the Surface class.
 */

#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "Device.h"
#include "Scalars.h"
#include "Vectors.h"
#include "DataIO.h"
#include "SHTrans.h"
#include "HelperFuns.h"
#include "Logger.h"

#include "SurfaceParams.h"
#include "OperatorsMats.h"

template <typename T, enum DeviceType DT> 
class Surface
{
  private:
    Device<DT> *device_;
    SurfaceParams<T> params_;
    
    SHTrans<T,DT> sht_;
    int max_n_surfs_;

  public:
    Vectors<T,DT> x_;
    Vectors<T,DT> normal_;
    Vectors<T,DT> cu_, cv_;
    Scalars<T,DT> w_;
    Scalars<T,DT> h_;
    Scalars<T,DT> k_;
    Vectors<T,DT> bending_force_;
    
    Scalars<T,DT> tension_;
    Vectors<T,DT> tensile_force_;

    T *all_rot_mats_;
    T *quad_weights_;

    Surface(Device<DT> *device_in, SurfaceParams<T> params_in, const OperatorsMats<T> &mats);
    ~Surface();

    void UpdateFirstForms();
    void UpdateAll();
    void Reparam();
    void SurfGrad(const Scalars<T,DT> &f_in, Vectors<T,DT> &grad_f_out);
    void SurfDiv(const Vectors<T,DT> &f_in, Scalars<T,DT> &div_f_out);
    void StokesMatVec(const Vectors<T,DT> &density_in, Vectors<T,DT> &velocity_out);
    void GetTension(const Vectors<T,DT> &v_in, const Vectors<T,DT> &v_ten_in, T *tension_out);
    T Area();
    void Volume();
    void Resize(int n_surfs_in);
    void Populate(const T *centers);
    bool IsAccurate();
    T* GetCenters(T* cnts);
    void UpdateNormal();

    //Work space
    Scalars<T,DT> S1, S2, S3, S4, S5, S6;
    Vectors<T,DT> V1, V2;
    T *shc, *work_arr, *alpha_p;
    T max_init_area_;
    Scalars<T,DT> S10;
    Vectors<T,DT> V10, V11, V12, V13;
    T *alpha_q;
    //Rotation
    Scalars<T,DT> w_sph_;
    T *rot_mat;
    T *sing_quad_weights_;
    double StokesMatVec_time_;
};

#include "Surface.cc"

#endif //_SURFACE_H_
