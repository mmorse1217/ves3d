/**
 * @file   vesicle.h
 * @author Rahimian, Abtin <abtin@romario>
 * @date   Mon Jan 18 14:08:45 2010
 * 
 * @brief  The header function of the class vesicle.
 *
 */

/**
 * @class vesicle
 * @brief The main class
 *
 * vesicle is the class used to ...
 * 
 */
template <typename T> class vesicle
{
public:
    int nv;/**< The number of vesicles */
    int p;/**< The spherical harmonics basis degree.*/
    int np;/**< number of discrization points on each vesicle */
    
    T* posVec;/**< The position vector for the grid points ordered as
	       * [x y z] first on each vesicle then across vesicles*/

    T** scaGeoProp;/**< The vector of scalar valued geometric properties, ordered as {W,H,K,...} */
    T** vecGeoProp/**< The vector of vector valued geometric  properties */

    int setPosVec(T* posIn);/**< Takes care of updating geoProps */

    vesicle();
    ~vesicle();

    int linCurv(T* vec, T* curv); 
    /**< Linear curvature operator */
    int Grad(T* vec, T* gradVec);
    /**< Surface gradient */
    int Div(T* vec, T* divVec);
    /**< surface divergence */
    int bendingOp(T* vec, T* bendForce);
    /**< Bending Operator */
    int tensionOp(T* vec, T* tensileForce);
    /**< Tension Operator */
};

int stokesMatVec(int p, int nv, T* posVec,T* W, T* density, T* sf); 
/// Self interaction stokes matVec

int DmSH(T* vecIn,int du,int dv, T* vecOut);
/// Spherical harmonic differentiation

int shAna(T* vecIn, T* vecOut);
/// Spherical harmonics transform

int shSyn(T* vecIn, T* vecOut);
/// Spherical harmonics inverse transform

int interpsh(T* vecIn, T* vecOut);
/// Interpolation function via spherical harmonics

int filtersh(T* vecIn, T* vecOut);
/// Filtering via spherical harmonics

int reparam(vecicle* ves);


