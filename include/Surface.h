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
#include <iostream>

template <typename T>
struct SurfaceParams 
{
    int p_;
    int n_surfs_;
    T kappa_;
    int filter_freq_;
    T rep_ts_;
    T rep_max_vel_;
    int rep_iter_max_;
    int rep_up_freq_;
    int rep_filter_freq_;
};

template<typename T>
ostream& operator<<(ostream& output, const SurfaceParams<T>& par)
{
    
    output<<"\n ------------------------------------"<<endl;
    output<<"  Surface properties"<<endl;
    output<<" ------------------------------------"<<endl;
    output<<"  p                 : "<<par.p_<<endl;
    output<<"  Number of surfaces: "<<par.n_surfs_<<endl;
    output<<"  kappa             : "<<par.kappa_<<endl;
    output<<"  filter_freq       : "<<par.filter_freq_<<endl;
    output<<"  rep_ts            : "<<par.rep_ts_<<endl;
    output<<"  rep_max_iter      : "<<par.rep_iter_max_<<endl;
    output<<"  rep_max_vel       : "<<par.rep_max_vel_<<endl;
    output<<"  rep_up_freq       : "<<par.rep_up_freq_<<endl;
    output<<"  rep_filter_freq   : "<<par.rep_filter_freq_<<endl;
    output<<" ------------------------------------"<<endl<<endl;
    
    return output;
}    

template <typename T> class Surface
{
  public:
    
    Device<T> &device_;

    
    SurfaceParams<T> params_;
    /** The vector holding the coordinates of the grid points,
     * supposing the Surface class, holds multiple surfaces, the order
     * will be \f$ \mathbf{x} = [X_1, Y_1, Z_1, \dots ,X_n, Y_n, Z_n],
     * \f$ where, for instance \f$X_1\f$ is x coordinate of the grid
     * points on the first surface; on each surface we sweep through
     * the points in a latitude major fashion.
     */
    Vectors<T> x_;

    /** 
     * The normal vector to the surface, that is given by \f[
     * \mathbf{n} = \frac{\mathbf{x_u \times x_v}}{W}, \f] where
     * \f$(u,v)\f$ are the parameters and W is the area element.
     */
    Vectors<T> normal_;

    
    Vectors<T> cu_, cv_;
    
    /** 
     * The area element size at the grid point, defined as \f[ W =
     * \sqrt{EG-F^2}, \f] where \f$E, F, G\f$ are the first
     * fundamental form coefficients.
     */
    Scalars<T> w_;

    /** 
     * The mean curvature at grid points, that is calculated according
     * to the formula \f[H = \frac{EN-2FM+GL}{2W^2},\f] where \f$E, F,
     * G\f$ are the first fundamental form coefficients and \f$L, M,
     * N\f$ are the second fundamental form coefficients.
     */
    Scalars<T> h_;

    /** 
     * The Gaussian curvature at grid points, that is calculated
     * according to the formula \f[K = \frac{LN-M^2}{W^2},\f] where
     * \f$E, F, G\f$ are the first fundamental form coefficients and
     * \f$L, M, N\f$ are the second fundamental form coefficients.
     */
    Scalars<T> k_;

    Vectors<T> bending_force_;

    T *tension_;
    Vectors<T> tensile_force_;
    
    //Rotation matrix
    T *all_rot_mats_;

    //Quadrature weights
    T* quad_weights_;
    
    Surface(Device<T> &device_in);
    Surface(Device<T> &device_in, SurfaceParams<T> params_in);
    Surface(Device<T> &device_in, SurfaceParams<T> params_in, const Vectors<T> &x_in);
    ~Surface();
    
    /** 
     * Setter method for the member function x_, although x_ is
     * public, setting it through this method envokes a call to
     * UpdateProps() that update the geometric properties of the
     * surface accordingly.
     */
    void SetX(const Vectors<T> &x_in);

    /// Recalculates the geometric properties of the surface.
    void UpdateFirstForms();
    void UpdateAll();

    void Reparam();

    /// The surface gradient operator.
    void SurfGrad(const Scalars<T> &f_in, Vectors<T> &grad_f_out);
    
    ///The surface divergence operator.
    void SurfDiv(const Vectors<T> &f_in, Scalars<T> &div_f_out);

    void StokesMatVec(const Vectors<T> &density_in, Vectors<T> &velocity_out);

    void GetTension(const Vectors<T> &v_in, const Vectors<T> &v_ten_in, T *tension_out);

    void Area();
    void Volume();

    void Resize(int n_surfs_in);

  public://private:
    //Work space
    Scalars<T> S1, S2, S3, S4, S5, S6;
    Vectors<T> V1, V2;
    T *shc, *work_arr, *alpha_p;
    
    //Work vectors for the up-sampling
    Vectors<T> V10, V11, V12, V13;
    Scalars<T> S10;
    T *alpha_q;

    void UpdateNormal();

    //Rotation
    Scalars<T> w_sph_;
    T *rot_mat;
    T *sing_quad_weights_;
    T *vel;
};



#include "Surface.cc"
#endif //_SURFACE_H_