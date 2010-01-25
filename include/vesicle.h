/// Work space manager
template <typename ScalarType> class WorkSpaceManager
{
  private:
  public:
    WorkSpaceManager();
    ~WorkSpaceManager();
    int PopWorkMemery();
    int PushWorkMemery();
};

template <typename ScalarType> SHScalars
{
 public:
	int number_of_scalar_Functions_;
	int p_;

    FieldsOnSphere();
    FieldsOnSphere(int num_of_fields_in, int num_colatitude_in, int num_longitude_in);
    FieldsOnSphere(int num_of_fields_in, int num_colatitude_in, int num_longitude_in, ScalarType *data_in);
    ~FieldsOnSphere();
    
    ScalarType* SetData(ScalarType *data_in, int data_length_in); 
    ScalarType* GetData();
};
	
template <typename ScalarType> SHVectors
{
private:
  public:
		FieldType X() data_ ;
		FieldType Y() data_ + length ; //etc..

    SetData();

		//gb    D3Vectors* MatVecProduct(MatType* A, D3Vectors *b); // remove
    
};
friend    CrossProduct(D3Vectors *aIn, D3Vectors *bIc, D3vectors *aCb);
//friend   operator+();  //gb avoid mallocs() problematic for large length
//friend   operator-();
//friend    operator=();
friend    DotProduct();

friend		axpy(a,x,y,c);       //gb ( c = a*x + y )
	
	


    

/// Spherical harmonics
template <typename FieldsOnSphere> class SphericalHarmonic
{  
 private:

    FieldsOnShpere 
 public:
    WorkSpaceManager *workVectorStack;
    int allDerivatives_FrequencyToPoints( VectorType *f, VectorType *Duf, *Dvf, *Duuf,...);
    int pointsToFrequency(VectorType *fin, VectorType *fout);
    int frequencyToPoints(VectorType *fin, VectorType *fout);
    int resample(VecType *fIn, VecType *fout);
    int resampleWithScaling(VecType *fIn, Real *Uscaling, Real *Vscaling, VecType *fout);
    
    // composite functions
    int allDerivatives_PointsToPoints( VectorType *f, VectorType *Duf, *Dvf, *Duuf,...);
};

/// General Genus Zero Surfaces
template <typename ScalarType> class Surface
{
private:
    FieldOnSphere<ScalarType> x_, h_, w_, k_, cu_, cv_, normal_;
    	
public:
    int p;
    int number_of_points;
    
    T* X(void);	       
    T* X(T *Xin);    
    
    Surface();
    ~Surface();

    int SurfGrad(T* vec, T* surf_grad_vec_out);
    int SurdDiv(T* vec, T* surf_div_out);
};

/// Vesicle specific
template <typename VectorType> class Vesicle : public Surface
{
    int ApproxCurvatureJacobian(VectorType* vec_in, VectorType* lininearized_curvature_out); 
    int BendingForce(VectorType* bending_force_out);
		int BendingForceMatVec(D3Vector in, D3Vector out);
    int TensileForce(VectorType* tension_in, VectorType* tensile_force_out);
};

int StokesMatVec(Vesicle *ves_in, T* surface_velocity_out);
int reparam(Vecicle* ves);

class StokesMatVec{
 puclic:
	init();
	matvec();
}

