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

template<typename SurfContainer, typename Interaction>
class InterfacialVelocity
{
  private:
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::device_type device_type;
    typedef typename SurfContainer::Sca_t Sca_t;
    typedef typename SurfContainer::Vec_t Vec_t;
    typedef OperatorsMats<Sca_t> Mats_t;

  public:
    InterfacialVelocity(SurfContainer &S_in, const Interaction &inter,
        Mats_t &mats, const Parameters<value_type> &params,
        const BgFlowBase<Vec_t> &bgFlow);
    ~InterfacialVelocity();

    Error_t updatePositionExplicit(const value_type &dt);
    Error_t updatePositionImplicit(const value_type &dt);
    Error_t reparam();

    Error_t getTension(const Vec_t &vel_in, Sca_t &tension) const;
    Error_t stokes(const Vec_t &force, Vec_t &vel) const;
    Error_t updateInteraction() const;

    Error_t operator()(const Vec_t &x_new, Vec_t &time_mat_vec) const;
    Error_t operator()(const Sca_t &tension, Sca_t &tension_mat_vec) const;

    void* usr_ptr_;
    Sca_t& tension(){ return tension_;}

  private:
    SurfContainer &S_;
    const Interaction &interaction_;
    const BgFlowBase<Vec_t> &bg_flow_;
    const Parameters<value_type> &params_;

    InterfacialForce<SurfContainer> Intfcl_force_;
    BiCGStab<Sca_t, InterfacialVelocity> linear_solver_;
    BiCGStab<Vec_t, InterfacialVelocity> linear_solver_vec_;

    value_type dt_;
    //Operators
    Sca_t w_sph_, w_sph_inv_;
    Sca_t sing_quad_weights_;
    Sca_t quad_weights_;
    Sca_t quad_weights_up_;

    SHTrans<Sca_t, SHTMats<value_type, device_type> > sht_;
    SHTrans<Sca_t, SHTMats<value_type, device_type> > sht_upsample_;

    mutable MovePole<Sca_t,Mats_t> move_pole;
    mutable Vec_t velocity_;
    mutable Sca_t tension_;

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
