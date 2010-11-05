#ifndef _BGFLOW_H_
#define _BGFLOW_H_

#include "BgFlowBase.h"
#define BGPI 3.1415926535897932

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
template<typename ScalarContainer, typename VecContainer>
class ParabolicFlow : public BgFlowBase<VecContainer>
{
  private:
    typedef typename VecContainer::value_type value_type;

  public:
    ParabolicFlow(value_type radius, value_type center_vel,
        value_type flow_dir_x = 1, value_type flow_dir_y = 0, 
        value_type flow_dir_z = 0);
    
    virtual void operator()(const VecContainer &pos, const value_type time, 
        VecContainer &vel_inf) const;
    
  private:
    void CheckContainers(const VecContainer &ref) const;

    value_type inv_radius2_, center_vel_;
    value_type flow_dir_x_, flow_dir_y_, flow_dir_z_;
    mutable VecContainer flow_direction_;
    mutable ScalarContainer s_wrk_;
};

////////////////////////////////////////////////////////////////////////////////
template<typename VecContainer>
class TaylorVortex : public BgFlowBase<VecContainer>
{
  private:
    typedef typename VecContainer::value_type value_type;

  public:
    explicit TaylorVortex(value_type strength = 1, value_type x_period = BGPI,
        value_type y_period = BGPI);
    
    virtual void operator()(const VecContainer &pos, const value_type time,
        VecContainer &vel_inf) const;
    
  private:
    value_type strength_, x_period_, y_period_;
    mutable VecContainer wrk_vec1_, wrk_vec2_;
};

#include "BgFlow.cc"

#endif
