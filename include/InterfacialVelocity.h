#ifndef _INTERFACIALVELOCITY_H_
#define _INTERFACIALVELOCITY_H_

#include "InterfacialForce.h"
#include "BiCGStab.h"
#include "Logger.h"
#include "Device.h"

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
    
  public:
    InterfacialVelocity(SurfContainer &S_in, Interaction &inter, 
        OperatorsMats<value_type, device_type> &mats, const Parameters<value_type> &params, 
        BackgroundFlow &bgFlow = ShearFlow<Vec_t>);
   
    void updatePositionExplicit(const value_type &dt);
    void updatePositionImplicit(const value_type &dt);
    void reparam();

    void getTension(const Vec_t &vel_in, Sca_t &tension) const;
    void stokes(const Vec_t &force, Vec_t &vel) const;
    void updateInteraction();

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
    
    //Workspace
    mutable Vec_t velocity, u1_, u2_, u3_;
    mutable Sca_t tension_, wrk_;

#ifndef NDEBUG
  public:
    bool benchmarkExplicit(Vec_t &Fb, Vec_t &SFb, Sca_t &tension, 
        Vec_t &vel, Vec_t &xnew, value_type tol) const;
        
    bool benchmarkBendingForce(const Vec_t &x, Vec_t &Fb, value_type tol) const;
    bool benchmarkStokes(const Vec_t &F, Vec_t &SF, value_type tol) const;
    bool benchmarkTension(const Vec_t &vel, Sca_t &tension, value_type tol) const;
    bool benchmarkTensileForce(const Sca_t &tension, Vec_t &Fs, value_type tol) const;
    
#endif //NDEBUG
};

#include "InterfacialVelocity.cc"

#endif // _INTERFACIALVELOCITY_H_
