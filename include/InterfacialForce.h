#ifndef _INTERFACIALFORCE_H_
#define _INTERFACIALFORCE_H_

template<typename ScalarContainer, typename VectorContainer,
         template<typename SC, typename VC> class SurfContainer>
class InterfacialForce
{
  public:
    void BendingForce(const SurfContainer<ScalarContainer, VectorContainer> &S, 
        VectorContainer &Fb) const;

    void TensileForce(const SurfContainer<ScalarContainer, VectorContainer> &S, 
        const ScalarContainer &tension, VectorContainer &Fs) const;
};

#include "InterfacialForce.cc"

#endif _INTERFACIALFORCE_H_
