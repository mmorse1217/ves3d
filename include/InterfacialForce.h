#ifndef _INTERFACIALFORCE_H_
#define _INTERFACIALFORCE_H_

#include "OperatorsMats.h"
#include "Parameters.h"
#include "SHTrans.h"
#include "SHTMats.h"
#include "Logger.h"
#include "HelperFuns.h"
#include "VesicleProps.h"

template<typename SurfContainer>
class InterfacialForce
{
  private:
    typedef typename SurfContainer::Sca_t Sca_t;
    typedef typename SurfContainer::Vec_t Vec_t;
    typedef typename SurfContainer::Arr_t Arr_t;
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::device_type device_type;
    typedef VesicleProperties<Arr_t> VProp_t;

    const Parameters<value_type> &params_;
    const VProp_t &ves_props_;
    SHTrans<Sca_t, SHTMats<value_type, device_type> > sht_   ;
    SHTrans<Sca_t, SHTMats<value_type, device_type> > sht_up_;

    mutable Vec_t cen;
    mutable Sca_t s1, s2;
    mutable Vec_t v1, ftmp;
    mutable SurfContainer* S_up;

  public:
    InterfacialForce(const Parameters<value_type> &params,
        const VProp_t &ves_porps,
        const OperatorsMats<Arr_t> &mats);
    ~InterfacialForce();

    void bendingForce(const SurfContainer &S, Vec_t &Fb) const;
    void linearBendingForce(const SurfContainer &S, const Vec_t &x, Vec_t &Fb) const;
    void linearBendingForcePerVesicle(const SurfContainer &S, const Vec_t &x, Vec_t &Fb, const int vesicle_i) const;
    void tensileForce(const SurfContainer &S, const Sca_t &tension, Vec_t &Fs) const;
    void gravityForce(const SurfContainer &S, const Vec_t &x, Vec_t &Fg) const;
    void gravityForcePerVesicle(const SurfContainer &S, const Vec_t &x, Vec_t &Fg, const int vesicle_i) const;

    void explicitTractionJump(
	const SurfContainer &S,
        Vec_t &F) const;

    void implicitTractionJump(
	const SurfContainer &S,
	const Vec_t &x,
	const Sca_t &tension,
        Vec_t &F) const;
    
    void implicitTractionJumpPerVesicle(
	const SurfContainer &S,
	const Vec_t &x,
	const Sca_t &tension,
        Vec_t &F, const int vesicle_i) const;
};

#include "InterfacialForce.cc"

#endif //_INTERFACIALFORCE_H_
