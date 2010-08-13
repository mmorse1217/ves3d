#ifndef _MONITOR_H_
#define _MONITOR_H_

#include "Logger.h"

template<typename EvolveSurface>
class MonitorBase{
  private:
    typedef typename EvolveSurface::value_type value_type;
        
  public:
    virtual ~MonitorBase();
    virtual bool operator()(const EvolveSurface *state, value_type &t, 
        value_type &dt) = 0;
};

//////////////////////////////////////////////////////////////////////////////////////////
template<typename EvolveSurface>
class Monitor : public MonitorBase<EvolveSurface>
{
  private:
    typedef typename EvolveSurface::value_type value_type;

    value_type time_hor_;
    bool save_flag_;
    value_type save_stride_;
    DataIO IO;
    mutable value_type A0, V0;
    
  public:
    Monitor(const Parameters<value_type> &params);
    ~Monitor();
    
    virtual bool operator()(const EvolveSurface *state, value_type &t, 
        value_type &dt);
};

#include "Monitor.cc"

#endif //_MONITOR_H_
