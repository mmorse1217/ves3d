#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "InterfacialVelocity.h"

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
    Monitor();
    bool operator()(const SurfContainer &state, 
        value_type &t, value_type &dt);
};

template<typename Container, typename Interaction>
class EvolveSurface
{
  private:
    typedef typename Container::value_type value_type;   
    OperatorsMats<value_type> &mats_;
    
  public:
    EvolveSurface(OperatorsMats<value_type> &mats);
    void operator()(Container &S, Interaction &Inter);
};

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
