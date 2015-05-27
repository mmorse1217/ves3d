#ifndef _INTERFACIALVELOCITY_H_
#define _INTERFACIALVELOCITY_H_

#include "Logger.h"
#include "Enums.h"
#include "Error.h"

#include "InterfacialForce.h"
#include "BiCGStab.h"
#include "SHTrans.h"
#include "Device.h"
#include <queue>
#include <memory>
#include "Enums.h"
#include "MovePole.h"
#include "BgFlowBase.h"
#include "OperatorsMats.h"
#include "ParallelLinSolverInterface.h"

#ifdef HAVE_PVFMM
#include "StokesVelocity.h"
#else
#include "StokesVelocitySeq.h"
#endif

template<typename SurfContainer, typename Interaction>
class InterfacialVelocity
{
  public:
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::device_type device_type;
    typedef typename SurfContainer::Arr_t Arr_t;
    typedef typename SurfContainer::Sca_t Sca_t;
    typedef typename SurfContainer::Vec_t Vec_t;
    typedef OperatorsMats<Arr_t> Mats_t;
    typedef InterfacialVelocity Matvec_t;
    typedef ParallelLinSolver<value_type> PSolver_t;
    typedef typename PSolver_t::matvec_type POp_t;
    typedef typename PSolver_t::vec_type PVec_t;
    typedef StokesVelocity<SurfContainer> Stokes_t;

    InterfacialVelocity(SurfContainer &S_in, const Interaction &inter,
        const Mats_t &mats, const Parameters<value_type> &params,
        const BgFlowBase<Vec_t> &bgFlow,
	PSolver_t *parallel_solver=NULL);

    ~InterfacialVelocity();

    Error_t Prepare(const SolverScheme &scheme) const;
    Error_t BgFlow(Vec_t &bg, const value_type &dt) const;

    Error_t AssembleRhsVel(PVec_t *rhs, const value_type &dt, const SolverScheme &scheme) const;
    Error_t AssembleRhsPos(PVec_t *rhs, const value_type &dt, const SolverScheme &scheme) const;
    Error_t AssembleInitial(PVec_t *u0, const value_type &dt, const SolverScheme &scheme) const;

    Error_t Solve(const PVec_t *rhs, PVec_t *u0, const value_type &dt, const SolverScheme &scheme) const;
    Error_t ConfigureSolver(const SolverScheme &scheme) const;
    Error_t ConfigurePrecond(const PrecondScheme &precond) const;
    Error_t Update(PVec_t *u0);

    Error_t updateJacobiExplicit   (const SurfContainer& S_, const value_type &dt, Vec_t& dx);
    Error_t updateJacobiGaussSeidel(const SurfContainer& S_, const value_type &dt, Vec_t& dx);
    Error_t updateJacobiImplcit    (const SurfContainer& S_, const value_type &dt, Vec_t& dx);
    Error_t updateImplicit         (const SurfContainer& S_, const value_type &dt, Vec_t& dx);

    Error_t reparam();

    Error_t getTension(const Vec_t &vel_in, Sca_t &tension) const;
    Error_t stokes(const Vec_t &force, Vec_t &vel) const;
    Error_t stokes_double_layer(const Vec_t &force, Vec_t &vel) const;
    Error_t updateFarField() const;

    Error_t EvaluateFarInteraction(const Vec_t &src, const Vec_t &fi, Vec_t &vel) const;
    Error_t CallInteraction(const Vec_t &src, const Vec_t &den, Vec_t &pot) const;

    Error_t operator()(const Vec_t &x_new, Vec_t &time_mat_vec) const;
    Error_t operator()(const Sca_t &tension, Sca_t &tension_mat_vec) const;

    Sca_t& tension(){ return tension_;}

  private:
    SurfContainer &S_;
    const Interaction &interaction_;
    const BgFlowBase<Vec_t> &bg_flow_;
    const Parameters<value_type> &params_;

    InterfacialForce<SurfContainer> Intfcl_force_;
    BiCGStab<Sca_t, InterfacialVelocity> linear_solver_;
    BiCGStab<Vec_t, InterfacialVelocity> linear_solver_vec_;

    // parallel solver
    PSolver_t *parallel_solver_;
    mutable bool psolver_configured_;
    mutable bool precond_configured_;
    mutable POp_t *parallel_matvec_;

    mutable PVec_t *parallel_rhs_;
    mutable PVec_t *parallel_u_;

    static Error_t ImplicitApply(const POp_t *o, const value_type *x, value_type *y);
    Error_t ImplicitMatvecPhysical(Vec_t &vox, Sca_t &ten) const;
    static Error_t ImplicitPrecond(const PSolver_t *ksp, const value_type *x, value_type *y);
    size_t stokesBlockSize() const;
    size_t tensionBlockSize() const;

    value_type dt_;

    Error_t EvalFarInter_Imp(const Vec_t &src, const Vec_t &fi, Vec_t &vel) const;
    Error_t EvalFarInter_ImpUpsample(const Vec_t &src, const Vec_t &fi, Vec_t &vel) const;

    //Operators
    Sca_t w_sph_, w_sph_inv_;
    Sca_t sing_quad_weights_;
    Sca_t quad_weights_;
    Sca_t quad_weights_up_;

    SHTrans<Sca_t, SHTMats<value_type, device_type> > sht_;
    SHTrans<Sca_t, SHTMats<value_type, device_type> > sht_upsample_;

    mutable Stokes_t stokes_;
    mutable MovePole<Sca_t,Mats_t> move_pole;
    mutable Vec_t pos_vel_;
    mutable Sca_t tension_;
    mutable Sca_t position_precond;
    mutable Sca_t tension_precond;
    mutable Arr_t vel_coeff_, dl_coeff_;

    //Workspace
    mutable std::queue<Sca_t*> scalar_work_q_;
    std::auto_ptr<Sca_t> checkoutSca() const;
    void recycle(std::auto_ptr<Sca_t> scp) const;
    mutable int checked_out_work_sca_;

    mutable std::queue<Vec_t*> vector_work_q_;
    std::auto_ptr<Vec_t> checkoutVec() const;
    void recycle(std::auto_ptr<Vec_t> vcp) const;
    mutable int checked_out_work_vec_;

    void purgeTheWorkSpace() const;

  //////////////////////////////////////////////////////////////////////////
  /// DEBUG mode methods ///////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
#ifndef NDEBUG
  public:
    bool benchmarkExplicit(Vec_t &Fb, Vec_t &SFb, Vec_t &vel,
        Sca_t &tension, Vec_t &xnew, value_type tol);

    bool benchmarkImplicit(Sca_t &tension, Vec_t &matvec,
        Vec_t &xnew, value_type tol);

    bool benchmarkBendingForce(const Vec_t &x, Vec_t &Fb, value_type tol) const;
    bool benchmarkStokes(const Vec_t &F, Vec_t &SF, value_type tol) const;
    bool benchmarkBgFlow(const Vec_t &SFb, Vec_t &vel, value_type tol) const;
    bool benchmarkTension(const Vec_t &vel, Sca_t &tension, value_type tol) const;
    bool benchmarkNewPostitionExplicit(Vec_t &xnew, value_type tol);

    bool benchmarkTensionImplicit(Sca_t &tension, value_type tol);
    bool benchmarkMatVecImplicit(const Vec_t &x, Vec_t &matvec, value_type tol);
    bool benchmarkNewPostitionImplicit(Vec_t &xnew, value_type tol);
#endif //NDEBUG
};

#include "InterfacialVelocity.cc"

#endif // _INTERFACIALVELOCITY_H_
