#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "TimeStepper.h"
#include "HelperFuns.h"
#include "InterfacialForce.h"
#include <iostream>
#include "DataIO.h"
#include "BiCGStab.h"

template<typename SurfContainer>
class InterfacialVelocity
{
  private:
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::Sca Sca;
    typedef typename SurfContainer::Vec Vec;
    
  public:
    InterfacialVelocity(SurfContainer *S_in);
    void operator()(const value_type &t, Vec &velocity) const;
    void operator()(const Sca &tension, Sca &div_stokes_fs) const;

    class TensionPrecond
    {
      public:
        void operator()(const Sca &in, Sca &out) const;
    } precond_;

  private:
    SurfContainer *S_;
    
    Sca w_sph_;
    Sca all_rot_mats_;
    Sca rot_mat_;
    Sca sing_quad_weights_;

    mutable Vec u1_, u2_;
    mutable Sca tension_;

    InterfacialForce<SurfContainer> Intfcl_force_;
    void GetTension(const Vec &vel_in, Sca &tension) const;
    void BgFlow(const Vec &pos, Vec &vel_Inf) const;
    void Stokes(const Vec &force, Vec &vel) const;
};

template<typename SurfContainer>
class Monitor
{
  private:
    typedef typename SurfContainer::value_type value_type;
    value_type time_hor_;
    
  public:
    Monitor(value_type th) : time_hor_(th) {};
    
    bool operator()(const SurfContainer &state, 
        value_type &t, value_type &dt)
    {
        typename SurfContainer::Sca A, V;
        A.replicate(state.getPosition());
        V.replicate(state.getPosition());
        state.area(A);
        state.volume(V);

        std::cout<<t<<" "<<A[0]<<"\t"<<V[0]<<std::endl;
        return(t<time_hor_);
    }       
};

template<typename SurfContainer, typename Forcing>
class ForwardEuler
{
  private:
    typedef typename SurfContainer::value_type value_type;
    typename SurfContainer::Vec velocity;
   public:
     void operator()(SurfContainer &S_in, value_type &t, value_type &dt, 
         Forcing &F, SurfContainer &S_out)
    {
        velocity.replicate(S_in.getPosition());
        F(t, velocity);
        axpy(dt, velocity, S_in.getPosition(), S_out.getPositionModifiable());

        ///@todo check for the norm of error as stopping criteria
        ///@todo upsampling?
        for(int ii=0; ii<10; ++ii)
        {
            S_out.getSmoothedShapePosition(velocity);
            axpy(static_cast<value_type>(-1), S_out.getPosition(), velocity, velocity);
            S_out.mapToTangentSpace(velocity);
            axpy(dt * 100, velocity, S_out.getPosition(), S_out.getPositionModifiable());
        }
    }
};

template<typename Container>
class EvolveSurface
{
  private:
    typedef typename Container::value_type value_type;
    
  public:
    
    void operator()(Container &S, value_type &time_horizon, value_type &dt)
    {
        value_type t(0);
        InterfacialVelocity<Container> F(&S);
        Monitor<Container> M(time_horizon);
        ForwardEuler<Container, InterfacialVelocity<Container> > U;
        
        TimeStepper<Container, 
            InterfacialVelocity<Container>, 
            ForwardEuler<Container, InterfacialVelocity<Container> >,
            Monitor<Container> > Ts;
    
        Ts(S, t, dt, F, U, M);
    }
};

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
