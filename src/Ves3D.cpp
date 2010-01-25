

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


function newTimeStep(){

	int nv; // number of vesicles;
	int  p; // discretization size
	int  m; // number of time steps
	double dt; 
	
	SHVectors vesPnts = InitVesiclesPos(n,p);
	Surface vesOld=Surface(vesPnts);
	Surface vesNew=Surface(ves0);
	SHVectors *force; // allocate memory

	sm = StokesMatVec();

	SHVectors work = SHScalar(vesPnts);
	
  for (i=0;i<m;i++){
		vesOld.BendingForce(force);
		sm.init(vesOld); 
		sm.matvec(force,work);
		axpy(dt,work,vesOld.X(),work);
		vesNew.X(work);
	}


		
		

	
 





