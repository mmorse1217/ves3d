#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "InterfacialVelocity.h"
#include "Logger.h"

template<typename SurfContainer>
class Monitor
{
  private:
    typedef typename SurfContainer::value_type value_type;
    value_type time_hor_;
    bool save_flag_;
    value_type save_stride_;
    DataIO<value_type, CPU> IO;
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
    OperatorsMats<value_type> &mats_;
    const Parameters<value_type> &params_;

  public:
    EvolveSurface(OperatorsMats<value_type> &mats, const Parameters<value_type> &params);
    void operator()(Container &S, Interaction &Inter);
};

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
