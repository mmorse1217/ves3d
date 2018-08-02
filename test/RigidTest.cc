#include <iostream>
#include "Vectors.h"
#include "Scalars.h"
#include "Surface.h"
#include "RigidParticle.h"
//#include "Parameters.h"
//#include "HelperFuns.h"
//#include <string>
#include "DataIO.h"
//#include "StokesVelocity.h"
#include "OperatorsMats.h"

const int DIM =3;
typedef Device<CPU> DevCPU;
extern const DevCPU the_cpu_device(0);

typedef double Real;
typedef Vectors<Real, DevCPU, the_cpu_device> Vec_t;
typedef Scalars<Real, DevCPU, the_cpu_device> Sca_t;
typedef Surface<Sca_t, Vec_t> Sur_t; 
typedef RigidParticle<Sur_t> RP;
typedef typename Sca_t::array_type Arr_t;
typedef OperatorsMats<Arr_t> Mats_t;

int main(int argc, char* argv[])
{
    VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
    DataIO myIO;
    int np, rank;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Comm_size(comm, &np);
    MPI_Comm_rank(comm, &rank);
    int nVec =1; // num particles
    
    Parameters<Real> sim_par;
    sim_par.sh_order = 8;
    int sh_order = sim_par.sh_order;
    sim_par.upsample_freq = 8;
    // initializing vesicle positions from text file
    Vec_t x0(nVec, sh_order);
    int fLen = x0.getStride();
    char fname[400];
    COUT("Loading initial shape");
    sprintf(fname, "precomputed/dumbbell_%d_double.txt",sim_par.sh_order);
    myIO.ReadData(FullPath(fname), x0, DataIO::ASCII, 0, fLen * DIM);

    COUT("Populating x0 (nves="<<nVec<<")");
    for(int ii=1;ii<nVec; ++ii)
        x0.getDevice().Memcpy(x0.getSubN_begin(ii),
                x0.begin(),
                x0.getTheDim() * fLen * sizeof(Real),
                Sur_t::device_type::MemcpyDeviceToDevice);

    //Reading operators from file
    COUT("Loading matrices");
    bool readFromFile = true;
    Mats_t mats(readFromFile, sim_par);

    //Creating objects
    COUT("Creating the surface object");
    Sur_t S(x0.getShOrder(),mats, &x0);

    COUT("Creating rigid particle");
    RP rp(S, sh_order, sh_order, sh_order);

    /*S.getPosition();
    Vec_t position(1, S.getPosition().getShOrder());
    f_i.getDevice().Memcpy(position.begin(), &(S_.getPosition().begin()),
                    vsz*sizeof(value_type), device_type::MemcpyDeviceToDevice);*/
    COUT("allocating vec");
    pvfmm::Vector<Real> position_pvfmm(rp.getSurface().getPosition().begin(),rp.getSurface().getPosition().size(), false);
    COUT("writign");

    rp.Solve();

    COUT("Density");
    for (int i = 0; i < rp.density_.size(); i++) {
        COUT(rp.density_.begin()[i]);
    }
    COUT("Translational velocity");
    for (int i = 0; i < rp.t_vel_.size(); i++) {
        COUT(rp.t_vel_.begin()[i]);
    }
    COUT("Rotational velocity");
    for (int i = 0; i < rp.r_vel_.size(); i++) {
        COUT(rp.r_vel_.begin()[i]);
    }
    //WriteVTK<Real>(S, 8, 8, "rigid_particle_test.vtp", 0, NULL, MPI_COMM_WORLD);

    WriteVTK<Sur_t>(S, "rigid_particle_test.vtp", MPI_COMM_WORLD, NULL, 32, 0);



    //for (int i = 0; i < 200; i++) {
    //rp.SetBoundaryData(?);
    //rp.Solve()
    //rp.update_position();
    //
    //}



    VES3D_FINALIZE();
    return 0;
}
