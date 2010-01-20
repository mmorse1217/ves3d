

function timeStepper(vesicle* ves, char* method, options* optsIn)
{
    
    T *bendForce;
    T *rhs1, *rhs2;
    T *sig;
    T *tensForce;
    int p, nv;
    double ts;
    
    ves.bendingOp(ves.scaGeoProp[1], bendForce);
    stokesMatVec(p, nv, S.posVec, S.scaGeoProp[2], bendForce, rhs1); 

    // Calculating tension
    GMRES(TENSIONMATVECPOINTER, rhs, sig);
    S.tensionOp(sig, tensForce);
    stokesMatVec(p, nv, S.posVec, S.scaGeoProp[2], tensForce, rhs2); 
    
    ves.posVec = ves.posVec + ts*(rhs1+rhs2); /// assuming overloading
					      /// the addition
    
}
