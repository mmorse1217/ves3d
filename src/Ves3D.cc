

// function timeStepper(vesicle* ves, char* method, options* optsIn)
// {
    
//     T *bendForce;
//     T *rhs1, *rhs2;
//     T *sig;
//     T *tensForce;
//     int p, nv;
//     double ts;

		
    
//     ves.bendingOp(ves.scaGeoProp[1], bendForce);
//     stokesMatVec(p, nv, S.posVec, S.scaGeoProp[2], bendForce, rhs1); 

//     // Calculating tension
//     GMRES(TENSIONMATVECPOINTER, rhs, sig);
//     S.tensionOp(sig, tensForce);
//     stokesMatVec(p, nv, S.posVec, S.scaGeoProp[2], tensForce, rhs2); 
    
//     ves.posVec = ves.posVec + ts*(rhs1+rhs2); /// assuming overloading
// 					      /// the addition
    
// }


// function newTimeStep(){

// 	int nv; // number of vesicles
// 	int  p; // discretization size
// 	int  m; // number of time steps
// 	double dt; 
	
//     int vec_size = 6*p*(p+1)*nv;
// 	SHVectors ves_pnts = InitVesiclesPos(nv,p);
// 	Vesicle ves(ves_pnts);
// 	SHVectors force(vec_size);
//     SHVectors work(vec_size), work2(vec_size);

//     sm = StokesMatVec();

//     for (i=0;i<m;i++){
// 		sm.init(ves); 
// 		ves.BendingForce(force);
// 		sm.matvec(force,work);
 
//         //Solve for tension;
//         GMRES(LinearizedTensionOperator(),*work, *ves.tension_); //Stokes should be included
//         TensileForce(force);
//         sm.matvec(force,work2);
//         work +=work2;
//         axpy(dt,work,vesOld.X(),work);
// 		ves.X(work);
// 	}
// }


		
		

	
 





