#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "Logger.h"
#include "Repartition.h"
#include "InterfacialVelocity.h"

template<typename SurfContainer>
class EvolveSurface
{
  private:
    typedef typename SurfContainer::value_type value_type;   
    typedef typename SurfContainer::device_type device_type;   
    typedef Monitor<SurfContainer> Mntr_t;
    typedef Interaction<value_type> Interaction_t;
    typedef RepartitionGateway<value_type> Repart_t;
    typedef InterfacialVelocity<SurfContainer, Interaction_t> IntVel_t;

    typedef void (IntVel_t::*Scheme_t)(const value_type &);

    OperatorsMats<typename Container::Sca_t> &mats_;
    const Parameters<value_type> &params_;
    Mntr &monitor_;
    const Repart &repartition_;

  public:
    //setters
    void setMonitor();
    void setInteraction();
    void setRepart();
    void setBgFlow();
    void setParameters();

    EvolveSurface(OperatorsMats<typename Container::Sca_t>
        &mats, const Parameters<value_type> &params, 
        Mntr &monitor, Repart &repartion);
    
    void operator()(Container &S, Interaction &Inter, 
        void* user_ptr = NULL);
};

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
