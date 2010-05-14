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

#include "OperatorsMats.h"
#include "HelperFuns.h"
#include "Logger.h"
#include <map>

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
    
    SurfaceParams();
    void SetMember(string var_name, string var_val);

  private:
    map<string, int> mapStringValues;

};

template <typename T, enum DeviceType DT> class Surface
{
  public:
    
    Device<DT> *device_;
    SurfaceParams<T> params_;
    SHTrans<T,DT> sht_;
    
    int max_n_surfs_;
    
    /** The vector holding the coordinates of the grid points,
     * supposing the Surface class, holds multiple surfaces, the order
     * will be \f$ \mathbf{x} = [X_1, Y_1, Z_1, \dots ,X_n, Y_n, Z_n],
     * \f$ where, for instance \f$X_1\f$ is x coordinate of the grid
     * points on the first surface; on each surface we sweep through
     * the points in a latitude major fashion.
     */
    Vectors<T,DT> x_;

    /** 
     * The normal vector to the surface, that is given by \f[
     * \mathbf{n} = \frac{\mathbf{x_u \times x_v}}{W}, \f] where
     * \f$(u,v)\f$ are the parameters and W is the area element.
     */
    Vectors<T,DT> normal_;

    
    Vectors<T,DT> cu_, cv_;
    
    /** 
     * The area element size at the grid point, defined as \f[ W =
     * \sqrt{EG-F^2}, \f] where \f$E, F, G\f$ are the first
     * fundamental form coefficients.
     */
    Scalars<T,DT> w_;

    /** 
     * The mean curvature at grid points, that is calculated according
     * to the formula \f[H = \frac{EN-2FM+GL}{2W^2},\f] where \f$E, F,
     * G\f$ are the first fundamental form coefficients and \f$L, M,
     * N\f$ are the second fundamental form coefficients.
     */
    Scalars<T,DT> h_;

    /** 
     * The Gaussian curvature at grid points, that is calculated
     * according to the formula \f[K = \frac{LN-M^2}{W^2},\f] where
     * \f$E, F, G\f$ are the first fundamental form coefficients and
     * \f$L, M, N\f$ are the second fundamental form coefficients.
     */
    Scalars<T,DT> k_;

    Vectors<T,DT> bending_force_;

    T *tension_;
    Vectors<T,DT> tensile_force_;
    
    //Rotation matrix
    T *all_rot_mats_;

    //Quadrature weights
    T *quad_weights_;
    
    Surface(Device<DT> *device_in, SurfaceParams<T> params_in, const OperatorsMats<T> &mats);
    ~Surface();
    
    /** 
     * Setter method for the member function x_, although x_ is
     * public, setting it through this method envokes a call to
     * UpdateProps() that update the geometric properties of the
     * surface accordingly.
     */
    void SetX(const Vectors<T,DT> &x_in);

    /// Recalculates the geometric properties of the surface.
    void UpdateFirstForms();
    void UpdateAll();

    void Reparam();

    /// The surface gradient operator.
    void SurfGrad(const Scalars<T,DT> &f_in, Vectors<T,DT> &grad_f_out);
    
    ///The surface divergence operator.
    void SurfDiv(const Vectors<T,DT> &f_in, Scalars<T,DT> &div_f_out);

    void StokesMatVec(const Vectors<T,DT> &density_in, Vectors<T,DT> &velocity_out);

    void GetTension(const Vectors<T,DT> &v_in, const Vectors<T,DT> &v_ten_in, T *tension_out);

    T Area();
    void Volume();

    void Resize(int n_surfs_in);

    void Populate(const T *centers);

    bool IsAccurate();

    T* GetCenters(T* cnts);

  public://private:
    //Work space
    Scalars<T,DT> S1, S2, S3, S4, S5, S6;
    Vectors<T,DT> V1, V2;
    T *shc, *work_arr, *alpha_p;
    
    T max_init_area_;
    //Work vectors for the up-sampling
    Scalars<T,DT> S10;
    Vectors<T,DT> V10, V11, V12, V13;

    T *alpha_q;

    void UpdateNormal();

    //Rotation
    Scalars<T,DT> w_sph_;
    T *rot_mat;
    T *sing_quad_weights_;

    double StokesMatVec_time_;
};



#include "Surface.cc"
#endif //_SURFACE_H_
