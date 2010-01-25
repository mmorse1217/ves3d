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

/// Cartesian vector
template <typename FieldType> class D3Vectors
{
  private:
    FieldType *X_,*Y_,*Z_;
    FieldType *data_;
  public:
    D3Vectors();
    ~D3Vectors();
    SetData();
    DotProduct();
    D3Vectors* CrossProduct(D3Vectors *a, D3Vectors *b);
    D3Vectors* MatVecProduct(MatType* A, D3Vectors *b);
    operator+();
    operator-();
    operator=();
};

// Spherical Calculus
/// Scalar fields defined on the unit sphere.
template <typename ScalarType> class FieldsOnSphere
{
private: 
    int number_of_fields_;
    int number_of_colatitude_grid_;
    int number_of_longitude_grid_;
    int fields_length_;
    ScalarType *data_;
    
public:
    FieldsOnSphere();
    FieldsOnSphere(int num_of_fields_in, int num_colatitude_in, int num_longitude_in);
    FieldsOnSphere(int num_of_fields_in, int num_colatitude_in, int num_longitude_in, ScalarType *data_in);
    ~FieldsOnSphere();
    
    ScalarType* SetData(ScalarType *data_in, int data_length_in); 
    ScalarType* GetData();
};

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

template <typename FieldType> class VectorsOnSphere : public D3Vectors
{
private:
    
  public:
    
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
    int BendingForce(VectorType* position_vector_in, VectorType* bending_force_out);
    int TensileForce(VectorType* tension_in, VectorType* tensile_force_out);
};

int StokesMatVec(Vesicle *ves_in, T* surface_velocity_out);
int reparam(Vecicle* ves);


