/**
 * @file   vesicle.h
 * @author Rahimian, Abtin <abtin@romario>
 * @date   Mon Jan 18 14:08:45 2010
 * 
 * @brief  The header function of the class vesicle.
 *
 */

/**
 * @class surface
 * @brief The main class
 *
 * surface is the class used to ...
 * 
 */

template <typename T> class vesicle{
    int bendingForce(T* vec, T* bendForce);
    /**< Bending Operator */
    int tensionForce(T* vec, T* tensileForce);
    /**< Tension Operator */
}



template <typename T> class surface
{
 private:
	VecType=.X_,E_,W_...
 
public:
    int nv;/**< The number of surfaces */
    int p;/**< The spherical harmonics basis degree.*/
    int np;/**< number of discrization points on each surface */
    
    T* X(void);/**< The position vector for the grid points ordered as
	       * [x y z] first on each vesicle then across vesicles*/
    T* X(T *Xin); // set calucaton and recalculate

		T* E(),W(),H,N,n,Cu,Cv;

    vesicle();
    ~vesicle();

    int approxCurvatureJacobian(T* vec, T* curv); 
    /**< Linearized curvature operator */
    int Grad(T* vec, T* gradVec);
    /**< Surface gradient */
    int Div(T* vec, T* divVec);
    /**< surface divergence */
};
int stokesMatVec(int p, int nv, T* posVec,T* W, T* density, T* sf); 
/// Self interaction stokes matVec

template <typename VectorType> class SH
{  
	
 public:
	VectorTypeManager *workVectorStack;
	int allDerivatives_FrequencyToPoints( VectorType *f, VectorType *Duf, *Dvf, *Duuf,...);
	int pointsToFrequency(VectorType *fin, VectorType *fout);
	int frequencyToPoints(VectorType *fin, VectorType *fout);
	int resample(VecType *fIn, VecType *fout);
	int resampleWithScaling(VecType *fIn, Real *Uscaling, Real *Vscaling, VecType *fout);

	// composite functions
	int allDerivatives_PointsToPoints( VectorType *f, VectorType *Duf, *Dvf, *Duuf,...);
};


template <typename ScalarType> class VecType
{
 public:
	int numberOfFunction;
	int order;  //p
	ScalarType *data;
	




int reparam(vecicle* ves);


