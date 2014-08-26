#ifndef _MONITOR_H_
#define _MONITOR_H_

#include "Logger.h"
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

    bool save_flag_;
    value_type save_stride_;
    DataIO IO;
    mutable value_type A0, V0;
    typename EvolveSurface::Sca_t area, vol;
    int last_save;
    const Parameters<value_type>& params_;

  public:
    Monitor(const Parameters<value_type> &params);
    ~Monitor();

    virtual Error_t operator()(const EvolveSurface *state, const value_type &t,
        value_type &dt);
};

#include "Monitor.cc"

#endif //_MONITOR_H_
