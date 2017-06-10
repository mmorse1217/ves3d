#ifndef _STOKES_VELOCITY_H_
#define _STOKES_VELOCITY_H_

#include <mpi.h>
#include "PVFMMInterface.h"
#include "NearSingular.h"
#include <matrix.hpp>

template <class Real>
class StokesVelocity{

    typedef typename pvfmm::Vector<Real> PVFMMVec;

  public:

    StokesVelocity(int sh_order, int sh_order_up, Real box_size=-1, Real repul_dist=1e-3, MPI_Comm c=MPI_COMM_WORLD);

    ~StokesVelocity();

    void SetSrcCoord(const PVFMMVec& S, int sh_order_up_self_=-1, int sh_order_up_=-1);

    template<class Vec>
    void SetSrcCoord(const Vec& S, int sh_order_up_self_=-1, int sh_order_up_=-1);

    void SetDensitySL(const PVFMMVec* force_single=NULL, bool add_repul=false);

    template<class Vec>
    void SetDensitySL(const Vec* force_single=NULL, bool add_repul=false);

    void SetDensityDL(const PVFMMVec* force_double=NULL);

    template<class Vec>
    void SetDensityDL(const Vec* force_double=NULL);

    void SetTrgCoord(const PVFMMVec* T);

    const PVFMMVec& operator()();
    const PVFMMVec& SelfInteraction();
    const PVFMMVec& FarInteraction();

    template<class Vec>
    void operator()(Vec& vel);
    template<class Vec>
    void SelfInteraction(Vec& vel);
    template<class Vec>
    void FarInteraction(Vec& vel);

    void setup_self();
    void setup_all();
    const Real* GetSLMatrixi(int i);
    const Real* GetDLMatrixi(int i);

    Real MonitorError(Real tol=1e-5);

    static void Test();

  private:

    StokesVelocity(const StokesVelocity &);
    StokesVelocity& operator=(const StokesVelocity &);

    int sh_order;
    int sh_order_up_self;
    int sh_order_up;
    Real box_size;
    MPI_Comm comm;


    PVFMMVec scoord;
    PVFMMVec scoord_far;
    PVFMMVec scoord_norm;
    PVFMMVec scoord_area;

    PVFMMVec force_single;
    PVFMMVec rforce_single; // force_single + repulsion
    PVFMMVec qforce_single; // upsample(rforce_single) * quadrature weights * area_element
    bool add_repul;

    PVFMMVec force_double;
    PVFMMVec uforce_double; // upsample(force_double)
    PVFMMVec qforce_double; // uforce_double * quadrature weights * area_element + normal

    bool trg_is_surf;
    PVFMMVec tcoord;
    PVFMMVec tcoord_repl;
    PVFMMVec trg_vel;


    // Self
    PVFMMVec SLMatrix, DLMatrix;
    PVFMMVec S_vel, S_vel_up;


    // Near
    NearSingular<Real> near_singular0; // Surface-to-Surface interaction
    NearSingular<Real> near_singular1; // Surface-to-Target interaction


    // Far
    bool fmm_setup;
    void* pvfmm_ctx;
    PVFMMVec fmm_vel;

};

#include "StokesVelocity.cc"

#endif // _STOKES_VELOCITY_H_
