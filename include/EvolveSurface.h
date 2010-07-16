#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "InterfacialVelocity.h"
#include "Logger.h"

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
        value_type &t, value_type &dt);
};

template<typename Container, typename Interaction>
class EvolveSurface
{
  private:
    typedef typename Container::value_type value_type;   
    typedef typename Container::device_type device_type;   
    typedef InterfacialVelocity<Container, Interaction> IntVel_t;
    typedef void (IntVel_t::*Scheme_t)(const value_type &);

    OperatorsMats<typename Container::Sca_t> &mats_;
    const Parameters<value_type> &params_;

  public:
    EvolveSurface(OperatorsMats<typename Container::Sca_t>
        &mats, const Parameters<value_type> &params);
    void operator()(Container &S, Interaction &Inter);
};

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
