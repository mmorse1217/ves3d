#ifndef _BGFLOW_H_
#define _BGFLOW_H_

#include "BgFlowBase.h"
#define M_PI 3.1415926535897932

/////////////////////////////////////////////////////////////////////////////////
template<typename VecContainer>
class ShearFlow : public BgFlowBase<VecContainer>
{
  private:
    typedef typename VecContainer::value_type value_type;

  public:
    ShearFlow(value_type shear_rate);  
    
    virtual void operator()(const VecContainer &pos, const value_type time,
        VecContainer &vel_inf) const;
    
  private:
    value_type shear_rate_;
};

////////////////////////////////////////////////////////////////////////////////
template<typename VecContainer>
class ParabolicFlow : public BgFlowBase<VecContainer>
{
  private:
    typedef typename VecContainer::value_type value_type;

  public:
    ParabolicFlow(value_type flow_curvature);
    
    virtual void operator()(const VecContainer &pos, const value_type time, 
        VecContainer &vel_inf) const;
    
  private:
    value_type flow_curvature_;
};

////////////////////////////////////////////////////////////////////////////////
template<typename VecContainer>
class TaylorVortex : public BgFlowBase<VecContainer>
{
  private:
    typedef typename VecContainer::value_type value_type;

  public:
    explicit TaylorVortex(value_type strength = 1, value_type x_period = 2 * M_PI,
        value_type y_period = 2 * M_PI);
    
    virtual void operator()(const VecContainer &pos, const value_type time,
        VecContainer &vel_inf) const;
    
  private:
    value_type strength_, x_period_, y_period_;
    mutable VecContainer wrk_vec1_, wrk_vec2_;
};

#include "BgFlow.cc"

#endif
