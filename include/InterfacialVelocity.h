#ifndef _INTERFACIALVELOCITY_H_
#define _INTERFACIALVELOCITY_H_

#include "Logger.h"
#include "Enums.h"
#include "Error.h"

#include "InterfacialForce.h"
#include "BiCGStab.h"
#include "BiCGStabTMP.h"
#include "SHTrans.h"
#include "Device.h"
#include <queue>
#include <memory>
#include "Enums.h"
#include "BgFlowBase.h"
#include "OperatorsMats.h"
#include "ParallelLinSolverInterface.h"
#include "VesicleProps.h"
#include "StokesVelocity.h"
#include "ContactInterface.h"
#include "GMRESLinSolver.h"
#include <unordered_map>
#include <map>
#include "VesBoundingBox.h"

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
    typedef StokesVelocity<value_type> Stokes_t;
    typedef SHTrans<Sca_t, SHTMats<value_type, device_type> > SHtrans_t;
    typedef VesicleProperties<Arr_t> VProp_t;

    InterfacialVelocity(SurfContainer &S_in, const Interaction &inter,
        const Mats_t &mats, const Parameters<value_type> &params,
        const VProp_t &ves_props, const BgFlowBase<Vec_t> &bgFlow,
        PSolver_t *parallel_solver=NULL);

    ~InterfacialVelocity();

    Error_t Prepare(const SolverScheme &scheme) const;
    Error_t BgFlow(Vec_t &bg, const value_type &dt) const;

    Error_t AssembleRhsVel(PVec_t *rhs, const value_type &dt, const SolverScheme &scheme) const;
    Error_t AssembleRhsPos(PVec_t *rhs, const value_type &dt, const SolverScheme &scheme) const;
    Error_t AssembleInitial(PVec_t *u0, const value_type &dt, const SolverScheme &scheme) const;
    Error_t ImplicitMatvecPhysical(Vec_t &vox, Sca_t &ten) const;

    Error_t Solve(const PVec_t *rhs, PVec_t *u0, const value_type &dt, const SolverScheme &scheme) const;
    Error_t ConfigureSolver(const SolverScheme &scheme) const;
    Error_t ConfigurePrecond(const PrecondScheme &precond) const;
    Error_t Update(PVec_t *u0);

    Error_t updateJacobiExplicit   (const SurfContainer& S_, const value_type &dt, Vec_t& dx);
    Error_t updateJacobiGaussSeidel(const SurfContainer& S_, const value_type &dt, Vec_t& dx);
    Error_t updateJacobiImplicit   (const SurfContainer& S_, const value_type &dt, Vec_t& dx);
    Error_t updateImplicit         (const SurfContainer& S_, const value_type &dt, Vec_t& dx);

    Error_t reparam();

    Error_t getTension(const Vec_t &vel_in, Sca_t &tension) const;
    Error_t stokes(const Vec_t &force, Vec_t &vel) const;
    Error_t stokesSLPerVesicle(const Vec_t &force, Vec_t &vel, const int vesicle_i) const;
    Error_t stokesDLPerVesicle(const Vec_t &force, Vec_t &vel, const int vesicle_i) const;
    Error_t stokes_double_layer(const Vec_t &force, Vec_t &vel) const;
    Error_t updateFarField() const;

    Error_t EvaluateFarInteraction(const Vec_t &src, const Vec_t &fi, Vec_t &vel) const;
    Error_t CallInteraction(const Vec_t &src, const Vec_t &den, Vec_t &pot) const;

    Error_t operator()(const Vec_t &x_new, Vec_t &time_mat_vec) const;
    Error_t operator()(const Sca_t &tension, Sca_t &tension_mat_vec) const;
    Error_t operator()(const Vec_t &x_new, const Sca_t &tension, 
        Vec_t &time_mat_vec, Sca_t &tension_mat_vec) const;
    Error_t operator()(const Vec_t &x_new, const Sca_t &tension, const Arr_t &lambda,
        Vec_t &time_mat_vec, Sca_t &tension_mat_vec, Arr_t &lambda_mat_vec) const;
    Error_t operator()(const Arr_t &lambda, Arr_t &lambda_mat_vec) const;
    Error_t operator()(const Vec_t &x_new, const Sca_t &tension, 
        Vec_t &time_mat_vec, Sca_t &tension_mat_vec, const int vesicle_i) const;

    
    Error_t CVJacobianTrans(const Arr_t &lambda, Vec_t &f_col) const;
    Error_t CVJacobian(const Vec_t &x_new, Arr_t &lambda_mat_vec) const;
    Error_t LCPSelect(const Arr_t &lambda, Arr_t &lambda_mat_vec) const;
    Error_t SolveLCP(Vec_t &u_lcp, Sca_t &ten_lcp, Arr_t &lambda_lcp, Arr_t &cvs) const;
    Error_t SolveLCPSmall(Arr_t &lambda_lcp, Arr_t &cvs) const;
    Error_t minmap(const Arr_t &xin1, const Arr_t &xin2, Arr_t &xout) const;
    Error_t projectU1(Vec_t &u1, const Vec_t &x_old) const;
    Error_t sca_abs(Sca_t &xin) const;
    Error_t TransferVesicle(std::vector<value_type> &pos_s, std::vector<value_type> &pos_e, 
            pvfmm::Vector<value_type> &pole_s_pole, pvfmm::Vector<value_type> &pole_e_pole) const;
    Error_t FormLCPMatrix(Arr_t &lcp_matrix) const;
    Error_t FormLCPMatrixSparse(Arr_t &lcp_matrix) const;
    Error_t GetDx(Vec_t &col_dx, Sca_t &col_tension, const Vec_t &col_f) const;
    Error_t GetColPos(const Vec_t &xin, std::vector<value_type> &pos_vec, pvfmm::Vector<value_type> &pos_pole) const;
    Error_t UpdateVgradInd(int *ind1, int *ind2, int base, size_t length) const;
    Error_t ParallelGetVolumeAndGradient(const Vec_t &X_s, const Vec_t &X_e) const;
    Error_t ParallelFormLCPMatrixSparse(std::map<std::pair<size_t, size_t>, value_type> &lcp_matrix) const;
    Error_t ParallelSolveLCPSmall(Arr_t &lambda, Arr_t &cvs) const;
    Error_t ParallelCVJacobianTrans(const Arr_t &lambda, Vec_t &f_col) const;
    Error_t ConfigureLCPSolver() const;
    Error_t ImplicitLCPMatvec(Arr_t &lambda) const;

    value_type StokesError(const Vec_t &x) const;

    Sca_t& tension(){ return tension_;}

    // contact
    mutable int num_cvs_;
    mutable std::vector<value_type> IV_;
    mutable Arr_t lcp_matrix_;
    mutable std::map<std::pair<size_t, size_t>, value_type> parallel_lcp_matrix_;
    mutable int current_vesicle_;

  private:
    SurfContainer &S_;
    const Interaction &interaction_;
    const BgFlowBase<Vec_t> &bg_flow_;
    const Parameters<value_type> &params_;
    const VProp_t &ves_props_;

    InterfacialForce<SurfContainer> Intfcl_force_;
    BiCGStab<Sca_t, InterfacialVelocity> linear_solver_;
    BiCGStab<Vec_t, InterfacialVelocity> linear_solver_vec_;
    BiCGStabTMP<Vec_t, Sca_t, InterfacialVelocity> linear_solver_vec_sca_;
    GMRESLinSolver<value_type> linear_solver_gmres_;

    // parallel solver
    PSolver_t *parallel_solver_;
    mutable bool psolver_configured_;
    mutable bool precond_configured_;
    mutable POp_t *parallel_matvec_;

    mutable PVec_t *parallel_rhs_;
    mutable PVec_t *parallel_u_;

    static Error_t ImplicitApply(const POp_t *o, const value_type *x, value_type *y);
    static Error_t ImplicitPrecond(const PSolver_t *ksp, const value_type *x, value_type *y);
    
    static Error_t JacobiImplicitApply(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y);
    static Error_t JacobiImplicitApplyPerVesicle(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y);
    static Error_t JacobiImplicitLCPApply(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y);
    static Error_t JacobiImplicitPrecond(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y);
    static Error_t JacobiImplicitPrecondPerVesicle(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y);
    static Error_t JacobiImplicitLCPPrecond(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y);
    static Error_t LCPApply(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y);
    static Error_t LCPPrecond(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y);
    static Error_t ImplicitLCPApply(const POp_t *o, const value_type *x, value_type *y);
    static Error_t ImplicitLCPPrecond(const PSolver_t *ksp, const value_type *x, value_type *y);
    
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

    SHtrans_t sht_;
    SHtrans_t sht_upsample_;

    mutable Stokes_t stokes_;
    mutable Vec_t pos_vel_;
    mutable Sca_t tension_;
    mutable Sca_t position_precond;
    mutable Sca_t tension_precond;

    //Workspace
    mutable SurfContainer* S_up_;
    mutable std::queue<Sca_t*> scalar_work_q_;
    std::auto_ptr<Sca_t> checkoutSca() const;
    void recycle(std::auto_ptr<Sca_t> scp) const;
    mutable int checked_out_work_sca_;

    mutable std::queue<Vec_t*> vector_work_q_;
    std::auto_ptr<Vec_t> checkoutVec() const;
    void recycle(std::auto_ptr<Vec_t> vcp) const;
    mutable int checked_out_work_vec_;

    void purgeTheWorkSpace() const;
    
    //Contact
    mutable ContactInterface CI_;
    mutable ContactInterface CI_pair_;
    mutable Vec_t vgrad_;
    mutable std::vector<int> vgrad_ind_;
    mutable std::map<int, Vec_t*> ghost_vgrad_;
    mutable std::map<int, std::vector<int>*> ghost_vgrad_ind_;
    mutable std::vector<int> PA_;
    mutable SurfContainer* S_i_;
    mutable std::vector<int> contact_vesicle_list_;
    mutable VesBoundingBox<value_type> *VBBI_;
};

#include "InterfacialVelocity.cc"

#endif // _INTERFACIALVELOCITY_H_
