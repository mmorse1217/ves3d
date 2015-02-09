#ifdef HAVE_PVFMM
#ifndef _STOKES_VELOCITY_H_
#define _STOKES_VELOCITY_H_

#include <mpi.h>
#include <vector.hpp>
#include <PVFMMInterface.h>

template<typename Surf_t>
class StokesVelocity{

    typedef typename Surf_t::Vec_t Vec_t;
    typedef typename Surf_t::value_type Real_t;
    typedef typename pvfmm::Vector<Real_t> PVFMMVec_t;

  public:

    StokesVelocity(MPI_Comm c=MPI_COMM_WORLD):comm(c),near_singular(c){
      pvfmm_ctx=PVFMMCreateContext<Real_t>();
      S=NULL;
      force_single=NULL;
      force_double=NULL;
      S_vel=NULL;
    }

    ~StokesVelocity(){
      PVFMMDestroyContext<Real_t>(&pvfmm_ctx);
    }

    void SetSrcCoord(const Surf_t& S_);
    void SetSurfaceVel(const Vec_t& S_vel_);
    void SetDensity(const Vec_t* force_single_=NULL, const Vec_t* force_double_=NULL);

    void SetTrgCoord(const Surf_t& T_);
    void SetTrgCoord(Real_t* trg_coord, size_t N);

    void operator()(Vec_t& T_vel);
    Real_t* operator()();

    static void Test();

  private:
    void* pvfmm_ctx;
    NearSingular<Surf_t> near_singular;

    const Surf_t* S;
    const Vec_t* force_single;
    const Vec_t* force_double;
    const Vec_t* S_vel;

    PVFMMVec_t T;
    PVFMMVec_t trg_vel;

    MPI_Comm comm;
};

#include "StokesVelocity.cc"

#endif // _STOKES_VELOCITY_H_
#endif
