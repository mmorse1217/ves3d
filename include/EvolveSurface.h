#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "TimeStepper.h"
#include "HelperFuns.h"
#include <iostream>

template<typename Container, typename SurfContainer>
class Stokes : public Forcing<Container>
{
  private:
    SurfContainer *S_;
    typename SurfContainer::Sca w_sph_;
    typename SurfContainer::Sca all_rot_mats_;
    typename SurfContainer::Sca rot_mat_;
    typename SurfContainer::Sca sing_quad_weights_;

  public:
    Stokes(SurfContainer *S_in);

    virtual void  operator()(const Container &interfacial_force, 
        const typename Container::value_type &t, 
        Container &velocity) const;
};

template<typename Container>
class SurfaceEvolMonitor : Monitor<Container>
{
  public:    
    virtual bool operator()(const Container &state, 
        const typename Container::value_type &t, 
        typename Container::value_type &dt) const
    {
        std::cout<<"Monitor"<<std::endl;
        return false;
    }
};

// template<typename Container, typename Forcing>
// class Discretization
// {
//   public:
//     typedef typename Container::value_type value_type;

//     virtual void operator()(const Container &old_state, 
//         const value_type &t, const Forcing &F, 
//         const value_type &dt, Container &new_state) const = 0;
// };

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
