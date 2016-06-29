#include <StokesVelocity.h>

int main(int argc, char** argv){
  MPI_Init(&argc, &argv);
  pvfmm::SetSigHandler();

  MPI_Comm comm=MPI_COMM_WORLD;
  pvfmm::Profile::Enable(true);

  typedef double Real;
  StokesVelocity<Real>::Test();

  pvfmm::Profile::print(&comm);
  MPI_Finalize();
}

