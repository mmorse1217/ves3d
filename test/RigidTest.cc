#include <iostream>
#include "Vectors.h"
#include "Scalars.h"
#include "Surface.h"
#include "RigidParticle.h"
//#include "Parameters.h"
//#include "HelperFuns.h"
//#include <string>
//#include "DataIO.h"
//#include "StokesVelocity.h"
//#include "OperatorsMats.h"

typedef Device<CPU> DevCPU;
extern const DevCPU the_cpu_device(0);

typedef double Real;
typedef Vectors<Real, DevCPU, the_cpu_device> Vec_t;
typedef Scalars<Real, DevCPU, the_cpu_device> Sca_t;
typedef Surface<Sca_t, Vec_t> Sur_t; 
typedef RigidParticle<Sur_t> RP;

int main(int argc, char* argv[])
{
  VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
  int np, rank;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &rank);

  RP rp(0,0,0);



  VES3D_FINALIZE();
  return 0;
}
