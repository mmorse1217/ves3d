#ifndef _INTERFACIALFORCE_H_
#define _INTERFACIALFORCE_H_

#include "Parameters.h"

template<typename SurfContainer>
class InterfacialForce
{
  private:
    typedef typename SurfContainer::Sca Sca;
    typedef typename SurfContainer::Vec Vec;

  public:
    void BendingForce(const SurfContainer &S, Vec &Fb) const;

    void TensileForce(const SurfContainer &S, const Sca &tension, Vec &Fs) const;
};

#include "InterfacialForce.cc"

#endif _INTERFACIALFORCE_H_
