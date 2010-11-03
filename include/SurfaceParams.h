#ifndef _SURFACEPARAMS_H_
#define _SURFACEPARAMS_H_

#include <map>

template <typename T>
struct SurfaceParams 
{
    int p_;
    int n_surfs_;
    T kappa_;
    int filter_freq_;
    T rep_ts_;
    T rep_max_vel_;
    int rep_iter_max_;
    int rep_up_freq_;
    int rep_filter_freq_;
    
    SurfaceParams();
    void SetMember(string var_name, string var_val);

  private:
    map<string, int> mapStringValues;

};

#include "SurfaceParams.cc"

#endif //_SURFACEPARAMS_H_
