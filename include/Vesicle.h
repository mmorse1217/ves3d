// /// Vesicle specific
// template <typename ScalarType> class Vesicle : public Surface
// {
//     public;
//     SHScalars tension_;
//     int ApproxCurvatureJacobian(SHVectors *vec_in, SHScalars *ininearized_curvature_out); 
//     int BendingForce(SHVectors *bending_force_out);
//     int TensileForce(SHVectors *tensile_force_out);

//     int LinearizedBendingOperator(SHScalars *curvature, SHVectors *bending_force_out);
//     int LinearizedTensionOperator(SHScalars *tension  , SHVectors *tensile_force_out);
//     //int BendingForceMatVec(SHVectors vec_in, SHVectors vec_out);
    
//     friend int StokesMatVec(Vesicle *ves_in,...); //arn?
// };

// class StokesMatVec{
//  puclic:
// 	init();
// 	matvec();
// };
