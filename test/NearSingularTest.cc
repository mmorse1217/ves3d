#include <StokesVelocity.h>

int main(int argc, char** argv){
  VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
  pvfmm::SetSigHandler();

  MPI_Comm comm=MPI_COMM_WORLD;
  pvfmm::Profile::Enable(true);

  typedef double Real;
  StokesVelocity<Real>::Test();

  pvfmm::Profile::print(&comm);
  VES3D_FINALIZE();

  return 0;
}
