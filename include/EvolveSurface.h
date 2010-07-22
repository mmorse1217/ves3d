#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "InterfacialVelocity.h"
#include "Logger.h"
#include "RepartionGateway.h"

template<typename SurfContainer>
class Monitor
{
  private:
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::Sca_t::device_type device_type;
    value_type time_hor_;
    bool save_flag_;
    value_type save_stride_;
    DataIO IO;
    value_type A0, V0;
    
  public:
    Monitor(const Parameters<value_type> &params);
    bool operator()(const SurfContainer &state, 
        value_type &t, value_type &dt, void* user = NULL);
};

template<typename Container, 
         typename Interaction, 
         typename Mntr = Monitor<Container>,
         typename Repart = RepartionGateway<typename Container::Sca_t> >
class EvolveSurface
{
  private:
    typedef typename Container::value_type value_type;   
    typedef typename Container::device_type device_type;   
    typedef InterfacialVelocity<Container, Interaction> IntVel_t;
    typedef void (IntVel_t::*Scheme_t)(const value_type &);

    OperatorsMats<typename Container::Sca_t> &mats_;
    const Parameters<value_type> &params_;
    Mntr monitor_;
    Repart repartion_;
  public:
    EvolveSurface(OperatorsMats<typename Container::Sca_t>
        &mats, const Parameters<value_type> &params, 
        Mntr monitor, Repart repartion = Repart() );
    
    void operator()(Container &S, Interaction &Inter);
};

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
