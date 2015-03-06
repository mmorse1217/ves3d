#ifndef _STOKES_VELOCITY_SEQ_H_
#define _STOKES_VELOCITY_SEQ_H_

#include "Parameters.h"
#include "OperatorsMats.h"
#include "MovePole.h"
#include <vector>

template<typename Surf_t>
class StokesVelocity{

    typedef typename Surf_t::Vec_t Vec_t;
    typedef typename Vec_t::scalars_type Sca_t;
    typedef typename Vec_t::array_type Arr_t;
    typedef OperatorsMats<Arr_t> Mats_t;
    typedef typename Surf_t::value_type Real_t;
    typedef typename Surf_t::device_type device_type;

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
        const Parameters<Real_t> &sim_par_);

    ~StokesVelocity();

    void SetSrcCoord(const Surf_t& S_);
    void SetDensitySL(const Vec_t* force_single_=NULL);
    void SetDensityDL(const Vec_t* force_double_=NULL);

    void SetTrgCoord(const Surf_t& T_);

    void operator()(Vec_t& T_vel, unsigned int flag=StokesVelocity::UpdateAll);
    Real_t* operator()(unsigned int flag=StokesVelocity::UpdateAll);

    static void Test();

  private:

    StokesVelocity(const StokesVelocity &);
    StokesVelocity& operator=(const StokesVelocity &);

    static void u_ref(const Real_t* coord, int n, Real_t* out);
    static void force(const Real_t* coord, int n, Real_t* out);

    const Surf_t* S;
    const Surf_t* T;
    const Vec_t* force_single;
    const Vec_t* force_double;

    Vec_t S_vel;
    Vec_t SL_vel;
    Vec_t DL_vel;
    std::vector<Real_t> trg_vel;

    int sh_order;
    Sca_t quad_weights_;
    Sca_t w_sph_, w_sph_inv_, sing_quad_weights_;

    MovePole<Sca_t,Mats_t> move_pole;
    const Parameters<Real_t> &sim_par;
};

#include "StokesVelocitySeq.cc"

#endif // _STOKES_VELOCITY_SEQ_H_
