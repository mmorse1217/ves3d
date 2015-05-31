#ifndef _STOKES_VELOCITY_H_
#define _STOKES_VELOCITY_H_

#include <mpi.h>
#include <vector.hpp>
#include "Parameters.h"
#include "OperatorsMats.h"
#include "PVFMMInterface.h"
#include "NearSingular.h"
#include "MovePole.h"

template<typename Surf_t>
class StokesVelocity{

    typedef typename Surf_t::Vec_t Vec_t;
    typedef typename Vec_t::scalars_type Sca_t;
    typedef typename Vec_t::array_type Arr_t;
    typedef OperatorsMats<Arr_t> Mats_t;
    typedef typename Surf_t::value_type Real_t;
    typedef typename Surf_t::device_type device_type;
    typedef typename pvfmm::Vector<Real_t> PVFMMVec_t;

  public:

    enum VelocityUpdateFlag{
      UpdateNone=0,
      UpdateSelf=1,
      UpdateNear=2,
      UpdateFar =4,
      UpdateAll =7
    };

    StokesVelocity(
        const OperatorsMats<Arr_t> &mats,
        const Parameters<Real_t> &sim_par_,
        MPI_Comm c=MPI_COMM_WORLD);

    ~StokesVelocity();

    void SetSrcCoord(const Surf_t& S_);
    void SetSurfaceVel(const Vec_t& S_vel_);
    void SetDensitySL(const Vec_t* force_single_=NULL);
    void SetDensityDL(const Vec_t* force_double_=NULL);

    void SetTrgCoord(const Surf_t& T_);
    void SetTrgCoord(Real_t* trg_coord, size_t N);

    void operator()(Vec_t& T_vel, unsigned int flag=StokesVelocity::UpdateAll);
    Real_t* operator()(unsigned int flag=StokesVelocity::UpdateAll);

    Real_t MonitorError();

    static void Test();

  private:

    StokesVelocity(const StokesVelocity &);
    StokesVelocity& operator=(const StokesVelocity &);

    void GetPole(const Vec_t& v,  PVFMMVec_t& pvfmm_v);
    void Vec2PVFMMVec(const Vec_t& v,  PVFMMVec_t& pvfmm_v);
    void Upsample(const Vec_t& v, Vec_t* v_out=NULL,  PVFMMVec_t* pvfmm_v=NULL);

    template<bool SL, bool DL>
    void SingularInteg();
    void Setup();

    const Vec_t& SelfInteraction(bool update_self);
    const PVFMMVec_t& NearInteraction(bool update_near);
    const PVFMMVec_t& FarInteraction(bool update_far);

    static void u_ref(const Real_t* coord, int n, Real_t* out);
    static void force(const Real_t* coord, int n, Real_t* out);


    const Surf_t* S;
    const Surf_t* S_far;
    Surf_t*       S_self;
    const Vec_t* force_single;
    const Vec_t* force_double;
    const Vec_t* S_vel_ptr;

    PVFMMVec_t src_coord_up;
    PVFMMVec_t qforce_single_up;
    PVFMMVec_t qforce_double_up;
    PVFMMVec_t force_double_up;
    PVFMMVec_t surf_vel_up;
    PVFMMVec_t trg_coord;
    bool trg_is_surf;

    Vec_t S_vel;
    Vec_t SL_vel;
    Vec_t DL_vel;
    PVFMMVec_t fmm_vel;
    PVFMMVec_t trg_vel;

    enum UpdateFlag{
      UpdateSrcCoord  = 1,
      UpdateTrgCoord  = 2,
      UpdateDensitySL = 4,
      UpdateDensityDL = 8,
      UpdateqDensitySL=16,
      UpdateqDensityDL=32,
      UpdateSurfaceVel=64,
    };
    unsigned int fmm_flag;
    unsigned int self_flag;

    int sh_order, sh_order_far, sh_order_self;
    Sca_t quad_weights_;              // for far-field integration (sh_order_far)
    Sca_t w_sph_inv_;                 // for singular integration  (sh_order)
    Sca_t sing_quad_weights_;         // for singular integration  (sh_order_self)
    Sca_t w_sph_sing_quad_weights_;   // for singular integration  (sh_order_self)
    PVFMMVec_t pole_quad;

    MPI_Comm comm;
    Real_t box_size_;
    MovePole<Sca_t,Mats_t> move_pole;
    const Parameters<Real_t> &sim_par;
    SHTrans<Sca_t, SHTMats<Real_t, device_type> > sht_   ;
    SHTrans<Sca_t, SHTMats<Real_t, device_type> > sht_up_;
    NearSingular<Surf_t> near_singular;
    void* pvfmm_ctx;
};

#include "StokesVelocity.cc"

#endif // _STOKES_VELOCITY_H_
