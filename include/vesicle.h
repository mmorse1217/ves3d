	
template <typename ScalarType> class FieldOnSphere
{
private: 
    int number_of_fields_;
    int each_field_length_;
    ScalarType *data_;

public:
    VectorsOnSphere();
    VectorsOnSphere(int num_of_fields_in, int each_field_length_in);
    VectorsOnSphere(int num_of_fields_in, int each_field_length_in, ScalarType *data_in);
    ~VectorsOnSphere();
    
    ScalarType* SetData(ScalarType *data_in,int data_length_in); 
}

template <typename ScalarType> class Surface
{
private:
    FieldOnSphere<ScalarType> x_, h_, w_, k_, cu_, cv_, normal_;
    	
public:
    int p;
    int number_of_points;
    
    T* X(void);	       
    T* X(T *Xin);    
    
    vesicle();
    ~vesicle();

    int SurfGrad(T* vec, T* surf_grad_vec_out);
    int SurdDiv(T* vec, T* surf_div_out);
};

template <typename VectorType> class Vesicle : public Surface<VectorType>{
    int ApproxCurvatureJacobian(VectorType* vec_in, VectorType* lininearized_curvature_out); 
    int BendingForce(VectorType* position_vector_in, VectorType* bending_force_out);
    int TensileForce(VectorType* tension_in, VectorType* tensile_force_out);
}


template <typename VectorsOnSphere> class SH
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

int StokesMatVec(Vesicle *ves_in, T* surface_velocity_out);
int reparam(Vecicle* ves);


