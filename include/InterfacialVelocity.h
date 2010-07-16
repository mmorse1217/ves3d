#ifndef _INTERFACIALVELOCITY_H_
#define _INTERFACIALVELOCITY_H_

#include "InterfacialForce.h"
#include "BiCGStab.h"
#include "Logger.h"
#include "Device.h"
#include <queue>
#include <memory>
#include "enums.h"
#include "MovePole.h"

template<typename VecContainer>
void ShearFlow(const VecContainer &pos, const typename 
    VecContainer::value_type shear_rate, VecContainer &vel_inf);


template<typename SurfContainer, 
         typename Interaction,
         typename BackgroundFlow = void(&)(const typename SurfContainer::Vec_t&, 
             typename SurfContainer::value_type, typename SurfContainer::Vec_t&) >
class InterfacialVelocity
{
  private:
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::device_type device_type;
    typedef typename SurfContainer::Sca_t Sca_t;
    typedef typename SurfContainer::Vec_t Vec_t;
    typedef OperatorsMats<Sca_t> Mats_t;
    
  public:
    InterfacialVelocity(SurfContainer &S_in, Interaction &inter, 
        Mats_t &mats, const Parameters<value_type> &params, 
        BackgroundFlow &bgFlow = ShearFlow<Vec_t>);
    ~InterfacialVelocity();
    
    void updatePositionExplicit(const value_type &dt);
    void updatePositionImplicit(const value_type &dt);
    void reparam();

    void getTension(const Vec_t &vel_in, Sca_t &tension) const;
    void stokes(const Vec_t &force, Vec_t &vel) const;
    void updateInteraction() const;
    
    void operator()(const Vec_t &x_new, Vec_t &time_mat_vec) const; 
    void operator()(const Sca_t &tension, Sca_t &tension_mat_vec) const;

  private:
    SurfContainer &S_;
    Interaction &interaction_;
    BackgroundFlow &bg_flow_;
    const Parameters<value_type> &params_;
    
    InterfacialForce<SurfContainer> Intfcl_force_;
    BiCGStab<Sca_t, InterfacialVelocity> linear_solver_;
    BiCGStab<Vec_t, InterfacialVelocity> linear_solver_vec_;

    value_type dt_;
    //Operators 
    Sca_t w_sph_;
    Sca_t all_rot_mats_;
    Sca_t rot_mat_;
    Sca_t sing_quad_weights_;
    Sca_t quad_weights_;

    mutable MovePole<Sca_t,Mats_t> move_pole;
    mutable Vec_t velocity;
    mutable Sca_t tension_;

    //Workspace
    mutable queue<Sca_t*> scalar_work_q_;
    auto_ptr<Sca_t> checkoutSca() const;
    void recycle(auto_ptr<Sca_t> scp) const;
    mutable int checked_out_work_sca_;
    
    mutable queue<Vec_t*> vector_work_q_;
    auto_ptr<Vec_t> checkoutVec() const;
    void recycle(auto_ptr<Vec_t> vcp) const;
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
