#ifndef _BGFLOW_H_
#define _BGFLOW_H_

#include "BgFlowBase.h"
#include "Parameters.h"

#define BGPI 3.1415926535897932

template<typename Vec_t>
Error_t BgFlowFactory(Parameters<typename Vec_t::value_type> &params, BgFlowBase<Vec_t> **vInf);

/////////////////////////////////////////////////////////////////////////////////
template<typename Vec_t>
class ShearFlowImp : public BgFlowBase<Vec_t>
{
  public:
    typedef typename Vec_t::value_type value_type;

    ShearFlowImp(value_type shear_rate);

    virtual void operator()(const Vec_t &pos, const value_type time,
        Vec_t &vel_inf) const;

  private:
    value_type shear_rate_;
};

////////////////////////////////////////////////////////////////////////////////
template<typename Vec_t>
class ParabolicFlowImp : public BgFlowBase<Vec_t>
{
  public:
    typedef typename Vec_t::value_type value_type;
    typedef typename Vec_t::scalars_type Sca_t;

    ParabolicFlowImp(value_type radius, value_type center_vel,
        value_type flow_dir_x = 1, value_type flow_dir_y = 0,
        value_type flow_dir_z = 0);

    virtual void operator()(const Vec_t &pos, const value_type time,
        Vec_t &vel_inf) const;

  private:
    void CheckContainers(const Vec_t &ref) const;

    value_type inv_radius2_, center_vel_;
    value_type flow_dir_x_, flow_dir_y_, flow_dir_z_;
    mutable Vec_t flow_direction_;
    mutable Sca_t s_wrk_;
};

////////////////////////////////////////////////////////////////////////////////
template<typename Vec_t>
class ExtensionalFlowImp : public BgFlowBase<Vec_t>
{
  public:
    typedef typename Vec_t::value_type value_type;
    typedef typename Vec_t::device_type DT;

    ExtensionalFlowImp(value_type rate);

    virtual void operator()(const Vec_t &pos,
        const value_type time,
        Vec_t &vel_inf) const;
  private:
    void AdjustCoeffs(int stride) const;

    value_type rate_;
    mutable Vec_t coeffs_;
};

////////////////////////////////////////////////////////////////////////////////
template<typename Vec_t>
class TaylorVortexImp : public BgFlowBase<Vec_t>
{
  public:
    typedef typename Vec_t::value_type value_type;

    explicit TaylorVortexImp(value_type strength = 1, value_type x_period = BGPI,
        value_type y_period = BGPI);

    virtual void operator()(const Vec_t &pos, const value_type time,
        Vec_t &vel_inf) const;

  private:
    value_type strength_, x_period_, y_period_;
    mutable Vec_t wrk_vec1_, wrk_vec2_;
};

////////////////////////////////////////////////////////////////////////////////
template<typename Vec_t>
class PeriodicFlowImp : public BgFlowBase<Vec_t>
{
  public:
    typedef typename Vec_t::value_type value_type;

    PeriodicFlowImp(value_type strength, value_type period=12);

    virtual void operator()(const Vec_t &pos, const value_type time,
        Vec_t &vel_inf) const;

  private:
    value_type period_;
    value_type strength_;
};

/////////////////////////////////////////////////////////////////////////////////
/*
 * only intended for testing reparametrization
 */
template<typename Vec_t>
class TwisterFlowImp : public BgFlowBase<Vec_t>
{
  public:
    typedef typename Vec_t::value_type value_type;

    TwisterFlowImp(value_type twist);

    virtual void operator()(const Vec_t &pos, const value_type time,
        Vec_t &vel_inf) const;

  private:
    value_type twist_;
};


#include "BgFlow.cc"
#endif
