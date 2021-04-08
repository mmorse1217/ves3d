#ifndef _MONITOR_H_
#define _MONITOR_H_

#include "Logger.h"
#include "Spharm.h"
#include "Enums.h"

template<typename EvolveSurface>
class MonitorBase{
  private:
    typedef typename EvolveSurface::value_type value_type;

  public:
    virtual ~MonitorBase();
    virtual Error_t operator()(const EvolveSurface *state, const value_type &t,
        value_type &dt) = 0;
};

//////////////////////////////////////////////////////////////////////////////////////////
template<typename EvolveSurface>
class Monitor : public MonitorBase<EvolveSurface>
{
  private:
    typedef typename EvolveSurface::value_type value_type;

    bool checkpoint_flag_;
    value_type checkpoint_stride_;
    DataIO IO_;
    //mutable value_type A0_, V0_;
    //typename EvolveSurface::Sca_t area0_, vol0_;
    typename EvolveSurface::Sca_t area_new, vol_new;
    int last_checkpoint_;
    int time_idx_;
    DictString_t d_;
    const Parameters<value_type> *params_;

  public:
    Monitor(const Parameters<value_type> *params);
    ~Monitor();
    mutable value_type A0_, V0_;
    typename EvolveSurface::Sca_t area0_, vol0_;

    virtual Error_t operator()(const EvolveSurface *state, const value_type &t,
        value_type &dt);
};

#include "Monitor.cc"

#endif //_MONITOR_H_
