#include <StokesVelocity.h>
#include <Vectors.h>

typedef Device<CPU> DevCPU;
extern const DevCPU the_cpu_device(0);

int main(int argc, char** argv){
  MPI_Init(&argc, &argv);
  pvfmm::SetSigHandler();

  MPI_Comm comm=MPI_COMM_WORLD;
  pvfmm::Profile::Enable(true);

  typedef double Real_t;
  typedef Vectors<Real_t, DevCPU, the_cpu_device> VecCPU_t;
  typedef typename VecCPU_t::scalars_type Sca_t;
  typedef Surface<Sca_t, VecCPU_t> Surf_t;
  StokesVelocity<Surf_t>::Test();

  pvfmm::Profile::print(&comm);
  MPI_Finalize();
}

