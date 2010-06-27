#ifndef _INTERFACIALVELOCITY_H_
#define _INTERFACIALVELOCITY_H_

#include "InterfacialForce.h"
#include "BiCGStab.h"
#include "Logger.h"

template<typename VecContainer>
void ShearFlow(const VecContainer &pos, const typename 
    VecContainer::value_type shear_rate, VecContainer &vel_inf);

template<typename SurfContainer, 
         typename Interaction,
         typename BackgroundFlow = void(&)(const typename SurfContainer::Vec&, 
             typename SurfContainer::value_type, typename SurfContainer::Vec&) >
class InterfacialVelocity
{
  private:
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::Sca Sca;
    typedef typename SurfContainer::Vec Vec;
    
  public:
    InterfacialVelocity(SurfContainer &S_in, Interaction &inter, 
        OperatorsMats<value_type, DataIO<value_type, CPU> > &mats, const Parameters<value_type> &params, 
        BackgroundFlow &bgFlow = ShearFlow<Vec>);
   
    void updatePositionExplicit(const value_type &dt);
    void updatePositionImplicit(const value_type &dt);
    void reparam();

    void operator()(const Vec &x_new, Vec &time_mat_vec) const; 
    void operator()(const Sca &tension, Sca &tension_mat_vec) const; 
  
  private:
    SurfContainer &S_;
    Interaction &interaction_;
    BackgroundFlow &bg_flow_;
    const Parameters<value_type> &params_;
    
    InterfacialForce<SurfContainer> Intfcl_force_;
    BiCGStab<Sca, InterfacialVelocity> linear_solver_;
    BiCGStab<Vec, InterfacialVelocity> linear_solver_vec_;

    value_type dt_;
    //Operators 
    Sca w_sph_;
    Sca all_rot_mats_;
    Sca rot_mat_;
    Sca sing_quad_weights_;
    Sca quad_weights_;
    
    //Workspace
    mutable Vec velocity, u1_, u2_, u3_;
    mutable Sca tension_, wrk_;

    void getTension(const Vec &vel_in, Sca &tension) const;
    void stokes(const Vec &force, Vec &vel) const;
    void updateInteraction();
};

#include "InterfacialVelocity.cc"

#endif // _INTERFACIALVELOCITY_H_
