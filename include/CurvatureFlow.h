#ifndef _CURVATUREFLOW_H_
#define _CURVATUREFLOW_H_

#include "TimeStepper.h"
#include <iostream>

template<typename SurfContainer>
class InterfacialVelocity
{
  private:
    typedef typename SurfContainer::value_type value_type;
    
  public:
    void  operator()(SurfContainer &S, value_type &t, 
        typename SurfContainer::Vec &velocity) const
    {
        xv(S.h_, S.normal_, velocity);
        axpy(.05, velocity, velocity);
    } 
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
        std::cout<<t<<std::endl;
        return(t<time_hor_);
    }       
};

template<typename SurfContainer, typename Forcing>
class ForwardEuler
{
  private:
    typedef typename SurfContainer::value_type value_type;
    
  public:
    void operator()(SurfContainer &S_in, value_type &t, value_type &dt, 
        Forcing &F, SurfContainer &S_out)
    {
        typename SurfContainer::Vec vel(S_in.x_.GetShOrder(), 
            S_in.x_.GetNumVecs());
        S_in.UpdateAll();
        F(S_in, t, vel);
        axpy(dt, vel, S_in.x_, S_out.x_);
    }
};

template<typename Container>
class CurvatureFlow
{
  private:
    typedef typename Container::value_type value_type;
  
public:

    void operator()(Container &S, value_type &time_horizon, value_type &dt)
    {
        value_type t = 0.0;
        InterfacialVelocity<Container> F;
        Monitor<Container> M(time_horizon);
        ForwardEuler<Container, InterfacialVelocity<Container> > U;
            
        TimeStepper<Container, 
            InterfacialVelocity<Container>, 
            ForwardEuler<Container, InterfacialVelocity<Container> >,
            Monitor<Container> > Ts;
        
        Ts(S, t, dt, F, U, M);
    }
};

#endif //_CURVATUREFLOW_H_
