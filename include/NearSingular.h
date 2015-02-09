#ifdef HAVE_PVFMM
#ifndef _NEAR_SINGULAR_H_
#define _NEAR_SINGULAR_H_

#include <mpi.h>
#include <vector.hpp>

template<typename Surf_t>
class NearSingular{

    typedef typename Surf_t::Vec_t Vec_t;
    typedef typename Surf_t::value_type Real_t;
    typedef typename pvfmm::Vector<Real_t> PVFMMVec_t;

  public:

    NearSingular(MPI_Comm c=MPI_COMM_WORLD):comm(c){};

    void SetSrcCoord(const Surf_t& S_);
    void SetSurfaceVel(const Vec_t& S_vel_);
    void SetDensity(const Vec_t* force_single_=NULL, const Vec_t* force_double_=NULL);

    void SetTrgCoord(const Surf_t& T_);
    void SetTrgCoord(const Real_t* trg_coord, size_t N);

    void EvaluateNearSingular(Vec_t& T_vel);
    Real_t* EvaluateNearSingular();

    static void TestNearSingular();

  private:

    void StokesNearSingular(Real_t r_near, const size_t* trg_cnt, const size_t* trg_dsp,  Real_t* trg_coord, Real_t* trg_veloc);

    const Surf_t* S;
    const Vec_t* force_single;
    const Vec_t* force_double;
    const Vec_t* S_vel;

    PVFMMVec_t T;
    PVFMMVec_t trg_vel;

    MPI_Comm comm;
};

#include "NearSingular.cc"

#endif // _NEAR_SINGULAR_H_
#endif
