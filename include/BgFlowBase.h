#ifndef _BGFLOWBASE_H_
#define _BGFLOWBASE_H_

/**
 * The abstract base class for the background flow. Any other
 * user-defined class for the backgound flow should inherit form this
 * class. As it is implicit in this class, the only assumption for the
 * background flow is the existence of the <tt>operator()</tt> that
 * given the Cartesian coordinate of the points and time, it returns
 * the velocity at those points.
 */
template<typename VecContainer>
class BgFlowBase
{
  public:
    /// The assumed interface.
    virtual void operator()(const VecContainer &position, const typename
	VecContainer::value_type time, VecContainer &v_inf) const = 0;
};
#endif
