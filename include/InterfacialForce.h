#ifndef _INTERFACIALFORCE_H_
#define _INTERFACIALFORCE_H_

#include "Parameters.h"
#include "Logger.h"

template<typename SurfContainer>
class InterfacialForce
{
  private:
    typedef typename SurfContainer::Sca_t Sca_t;
    typedef typename SurfContainer::Vec_t Vec_t;
    typedef typename SurfContainer::value_type value_type;

    value_type bending_modulus_;

    mutable Sca_t s1, s2;
    mutable Vec_t v1;

  public:
    InterfacialForce(const Parameters<value_type> &params);
    ~InterfacialForce();

    void bendingForce(const SurfContainer &S, Vec_t &Fb) const;
    void linearBendingForce(const SurfContainer &S, 
        const Vec_t &x_new, Vec_t &Fb) const;
    void tensileForce(const SurfContainer &S, const Sca_t &tension, 
        Vec_t &Fs) const;
};

#include "InterfacialForce.cc"

#endif //_INTERFACIALFORCE_H_
