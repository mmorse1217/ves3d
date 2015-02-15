#ifdef HAVE_PVFMM
#ifndef _STOKES_VELOCITY_H_
#define _STOKES_VELOCITY_H_

#include <mpi.h>
#include <vector.hpp>
#include <profile.hpp>

#include "Parameters.h"
#include "OperatorsMats.h"
#include "PVFMMInterface.h"
#include "NearSingular.h"

#include "DataIO.h"
#include "BgFlow.h"
#include "Array.h"
#include "Vectors.h"
#include "Surface.h"

#include "MovePole.h"
#include "VesInteraction.h"
#include "InterfacialVelocity.h"

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
        OperatorsMats<Arr_t> &mats,
        const Parameters<Real_t> &sim_par_,
        //const BgFlowBase<Vec_t> &bgFlow_,
        MPI_Comm c=MPI_COMM_WORLD):
      move_pole(mats),
      sim_par(sim_par_),
      //bg_flow(bgFlow_),
      near_singular(c),
      comm(c)
  {
      pvfmm_ctx=PVFMMCreateContext<Real_t>();
      S=NULL;
      force_single=NULL;
      force_double=NULL;
      S_vel_ptr=NULL;
      fmm_flag =StokesVelocity::UpdateNone;
      self_flag=StokesVelocity::UpdateNone;

      sh_order=sim_par.sh_order;
      { // quadrature weights
        quad_weights_.resize(1,sh_order);
        quad_weights_.getDevice().Memcpy(quad_weights_.begin(),
            mats.quad_weights_,
            quad_weights_.size() * sizeof(Real_t),
            Surf_t::device_type::MemcpyDeviceToDevice);
      }

      { // Set w_sph_, w_sph_inv_, sing_quad_weights_
        w_sph_.resize(1, sh_order);
        w_sph_inv_.resize(1, sh_order);
        int np = w_sph_.getStride();

        w_sph_.getDevice().Memcpy(w_sph_.begin(), mats.w_sph_,
            np * sizeof(Real_t), device_type::MemcpyDeviceToDevice);
        xInv(w_sph_,w_sph_inv_);

        //Singular quadrature weights
        sing_quad_weights_.resize(1,sh_order);
        sing_quad_weights_.getDevice().Memcpy(sing_quad_weights_.begin(),
            mats.sing_quad_weights_, sing_quad_weights_.size() *
            sizeof(Real_t),
            device_type::MemcpyDeviceToDevice);
      }
    }

    ~StokesVelocity(){
      PVFMMDestroyContext<Real_t>(&pvfmm_ctx);
    }

    void SetSrcCoord(const Surf_t& S_);
    void SetSurfaceVel(const Vec_t& S_vel_);
    void SetDensitySL(const Vec_t* force_single_=NULL);
    void SetDensityDL(const Vec_t* force_double_=NULL);

    void SetTrgCoord(const Surf_t& T_);
    void SetTrgCoord(Real_t* trg_coord, size_t N);

    void operator()(Vec_t& T_vel, unsigned int flag=StokesVelocity::UpdateAll);
    Real_t* operator()(unsigned int flag=StokesVelocity::UpdateAll);

    static void Test();

  private:

    StokesVelocity(const StokesVelocity &);
    StokesVelocity& operator=(const StokesVelocity &);

    void* pvfmm_ctx;
    NearSingular<Surf_t> near_singular;

    const Surf_t* S;
    const Vec_t* force_single;
    const Vec_t* force_double;
    const Vec_t* S_vel_ptr;

    PVFMMVec_t qforce_single;
    PVFMMVec_t qforce_double;
    Vec_t S_vel;

    PVFMMVec_t src_coord;
    PVFMMVec_t trg_coord;

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

    int sh_order;
    Sca_t quad_weights_;
    Sca_t w_sph_, w_sph_inv_, sing_quad_weights_;

    MovePole<Sca_t,Mats_t> move_pole;
    const Parameters<Real_t> &sim_par;
    //const BgFlowBase<Vec_t> &bg_flow;

    MPI_Comm comm;
};

#include "StokesVelocity.cc"

#endif // _STOKES_VELOCITY_H_
#endif
