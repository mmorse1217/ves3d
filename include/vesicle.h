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
 private:
    ScalarType *data_;
 public:
	int number_of_scalar_Functions_;
	int p_;

    SHScalars();
    ~SHScalars();
    
    ScalarType* SetData(ScalarType *data_in, int data_length_in); 
    ScalarType* GetData(ScalarType *data_out, int data_length);
};
	
template <typename ScalarType> SHVectors : public SHScalars
{
 private:
 public:
    friend SHVectors* operator+=(SHVectors *rhs);
    //friend   operator+(); //gb avoid mallows() problematic for large length
    //friend   operator-();
    //friend    operator=();
    friend CrossProduct(SHVectors *a_in, SHVectors *b_in, SHVectors *aCb_out);
    friend DotProduct(SHVectors *a_in, SHVectors *b_in, SHScalars *aDb_out);
    friend axpy(a,x,y,c); //gb ( c = a*x + y )    
};

/// Spherical harmonics
template <typename SHScalars> class SphericalHarmonic
{  
 private:
 public:
    WorkSpaceManager *workVectorStack;
    int AllDerivatives_FrequencyToPoints( SHScalars *f_in, SHScalars *Duf_out, SHScalars *Dvf_out, SHScalars *Duuf_out, SHScalars *Duvf_out, SHScalars *Dvvf_out);
    int PointsToFrequency(SHScalars *f_in, SHScalars *f_out);
    int FrequencyToPoints(SHScalars *f_in, SHScalars *f_out);
    int Resample(SHScalars *f_in, SHScalars *f_out);
    int ResampleWithScaling(SHScalars *f_in, int *p_scaling_in, SHScalars *f_out);
    
    // composite functions
};

/// General Surfaces
template <typename ScalarType> class Surface
{
private:
    SHVectors<ScalarType> x_, normal_;
    SHScalars<ScalarType> h_, w_, k_, cu_, cv_;
    
public:
    int p;
    
    //T* X(void);	       arn??
    //T* X(T *Xin);        arn??
    
    Surface();
    ~Surface();

    int SurfGrad(SHScalars *scalar_in, SHVectors* surf_grad_vec_out);
    int SurdDiv(SHVectors* vector_in, SHScalars* surf_div_out);
};

/// Vesicle specific
template <typename ScalarType> class Vesicle : public Surface
{
    public;
    SHScalars tension_;
    int ApproxCurvatureJacobian(SHVectors *vec_in, SHScalars *ininearized_curvature_out); 
    int BendingForce(SHVectors *bending_force_out);
    int TensileForce(SHVectors *tensile_force_out);

    int LinearizedBendingOperator(SHScalars *curvature, SHVectors *bending_force_out);
    int LinearizedTensionOperator(SHScalars *tension  , SHVectors *tensile_force_out);
    //int BendingForceMatVec(SHVectors vec_in, SHVectors vec_out);
    
    friend int StokesMatVec(Vesicle *ves_in,...); //arn?
};

class StokesMatVec{
 puclic:
	init();
	matvec();
}

