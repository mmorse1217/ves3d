#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "TimeStepper.h"
#include "HelperFuns.h"
#include "InterfacialForce.h"
#include <iostream>
#include "BiCGStab.h"
#include "Parameters.h"
#include "VesInteraction.h"
#include "GLIntegrator.h"
#include "SHTrans.h"
#include "OperatorsMats.h"

template<typename SurfContainer, typename Interaction>
class InterfacialVelocity
{
  private:
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::Sca Sca;
    typedef typename SurfContainer::Vec Vec;
    
  public:
    InterfacialVelocity(SurfContainer &S_in, Interaction &inter, 
        OperatorsMats<value_type> &mats);
    
    void operator()(const value_type &t, Vec &velocity) const;
    void operator()(const Sca &tension, Sca &div_stokes_fs) const;

    class TensionPrecond
    {
      public:
        void operator()(const Sca &in, Sca &out) const;
    } precond_;

    void reparam();

  private:
    SurfContainer &S_;
    
    Sca w_sph_;
    Sca all_rot_mats_;
    Sca rot_mat_;
    Sca sing_quad_weights_;
    Sca quad_weights_;
    
    mutable Vec u1_, u2_, u3_;
    mutable Sca tension_, wrk_;

    InterfacialForce<SurfContainer> Intfcl_force_;

    Interaction &interaction_;

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
    bool save_flag_;
    int save_stride_;
    DataIO<value_type, CPU> IO;
    
  public:
    Monitor() : 
        time_hor_(Parameters<value_type>::getInstance().time_horizon),
        save_flag_(Parameters<value_type>::getInstance().save_data),
        save_stride_(Parameters<value_type>::getInstance().save_stride),
        IO(SurfContainer::Vec::getDevice())
    {};
    
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

        F.reparam();
    }
};

template<typename Container, typename Interaction>
class EvolveSurface
{
  private:
    typedef typename Container::value_type value_type;   
    OperatorsMats<value_type> &mats_;
    
  public:
    EvolveSurface(OperatorsMats<value_type> &mats) :
        mats_(mats){}

    void operator()(Container &S, Interaction &Inter)
    {
        value_type t(0);
        value_type dt(Parameters<value_type>::getInstance().ts);
        
        InterfacialVelocity<Container, Interaction> F(S, Inter, mats_);
        
        Monitor<Container> M;
        ForwardEuler<Container, InterfacialVelocity<Container, Interaction> > U;
        
        TimeStepper<Container, 
            InterfacialVelocity<Container, Interaction>, 
            ForwardEuler<Container, InterfacialVelocity<Container, Interaction> >,
            Monitor<Container> > Ts;
    
        Ts(S, t, dt, F, U, M);
    }
};

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
