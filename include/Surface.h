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

template <typename T> class Surface
{
  public:
    
    Device<T> &device_;

    /// The spherical harmonics expansion truncation frequency.
    int p_;
    
    /// The number of separate surfaces in the Surface class.
    int n_surfs_;

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
        
    Surface(Device<T> &device_in);
    Surface(Device<T> &device_in, int p_in, int n_surfs_in);
    Surface(Device<T> &device_in, int p_in, int n_surfs_in, const Vectors<T> &x_in);
    
    /** 
     * Setter method for the member function x_, although x_ is
     * public, setting it through this method envokes a call to
     * UpdateProps() that update the geometric properties of the
     * surface accordingly.
     */
    void SetX(const Vectors<T> &x_in);

    /// Recalculates the geometric properties of the surface.
    void UpdateProps();

    /// The surface gradient operator.
    void SurfGrad(const Scalars<T> &f_in, Vectors<T> &grad_f_out);
    
    ///The surface divergence operator.
    void SurfDiv(const Vectors<T> &f_in, Scalars<T> &div_f_out);

    void StokesMatVec(const Vectors<T> &density_in, Vectors<T> &velocity_out);

    //Work space -- TO BE REMOVED
    Vectors<T> xu, xv, xuv, xvv, xuu, bending_force;
    Scalars<T> E, F, G, L, M, N, W2, temp;
};

#include "Surface.cc"
#endif //_SURFACE_H_
