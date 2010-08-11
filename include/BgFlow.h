#ifndef _BGFLOW_H_
#define _BGFLOW_H_

#define M_PI 3.1415926535897932

/**
 * The abstract base class for the background flow. Any other
 * user-defined class for the backgound flow should inherit form this
 * class. As it is implicit in this class, the only assumption for the
 * background flow is the existence of the <tt>operator()</tt> that
 * given the cartesian coordinate of the points, and time, it returns
 * the velocity at those points.
 */
template<typename VecContainer>
class BgFlow
{
  public:
    /// The assumed interface.
    virtual void operator()(const VecContainer &position, const typename 
        VecContainer::value_type time, VecContainer &v_inf) = 0;
};

template<typename VecContainer>
class ShearFlow : public BgFlow<VecContainer>
{
  private:
    typedef typename VecContainer::value_type value_type;

  public:
    ShearFlow(value_type shear_rate);  
    
    virtual void operator()(const VecContainer &pos, const value_type time,
        VecContainer &vel_inf) ;
    
  private:
    value_type shear_rate_;
};

template<typename VecContainer>
class ParabolicFlow : public BgFlow<VecContainer>
{
  private:
    typedef typename VecContainer::value_type value_type;

  public:
    ParabolicFlow(value_type flow_curvature);
    
    virtual void operator()(const VecContainer &pos, const value_type time, 
        VecContainer &vel_inf);
    
  private:
    value_type flow_curvature_;
};

template<typename VecContainer>
class TaylorVortex : public BgFlow<VecContainer>
{
  private:
    typedef typename VecContainer::value_type value_type;

  public:
    explicit TaylorVortex(value_type strength = 1, value_type x_period = 2 * M_PI,
        value_type y_period = 2 * M_PI);
    
    virtual void operator()(const VecContainer &pos, const value_type time,
        VecContainer &vel_inf);
    
  private:
    value_type strength_, x_period_, y_period_;
    VecContainer wrk_vec1_, wrk_vec2_;
};

#include "BgFlow.cc"

#endif
