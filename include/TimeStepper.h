#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

#include "Surface.h"

template<typename T> 
class TimeStepper
{
  public:
    T ts_;
    int n_steps_;
    Surface<T> &vesicle_;
    Vectors<T> velocity_, vel_tension_;
    

    TimeStepper(T ts_in, int n_steps_in, Surface<T> &ves_in) : 
        ts_(ts_in),
        n_steps_(n_steps_in), 
        vesicle_(ves_in),
        velocity_(vesicle_.device_, vesicle_.p_, vesicle_.n_surfs_),
        vel_tension_(vesicle_.device_, vesicle_.p_, vesicle_.n_surfs_)
    {};

    ~TimeStepper(){};

    void EvolveInTime()
    {
#ifndef NDEBUG
        cout<<" ------------------------------------"<<endl;
        cout<<"  Number of vesicles: "<<vesicle_.n_surfs_<<endl;
        cout<<"  p                 : "<<vesicle_.p_<<endl;
        cout<<"  ts                : "<<ts_<<endl;
        cout<<"  n_steps           : "<<n_steps_<<endl<<endl;

        cout<<"  kappa             : "<<vesicle_.kappa_<<endl;
        cout<<"  rep ts            : "<<vesicle_.rep_ts_<<endl;
        cout<<"  rep max iter      : "<<vesicle_.iter_max_<<endl;
        cout<<"  rep max vel       : "<<vesicle_.max_vel_<<endl;
        cout<<" ------------------------------------"<<endl<<endl;
#endif

        for(int idx=0;idx<n_steps_;idx++)
        {
            vesicle_.UpdateAll();
            //Interaction
            vesicle_.StokesMatVec(vesicle_.bending_force_, velocity_);
            vesicle_.StokesMatVec(vesicle_.tensile_force_, vel_tension_);
            //Calculate tension
            axpy((T) 1.0, velocity_, vel_tension_, velocity_);
            axpy(ts_, velocity_, vesicle_.x_, vesicle_.x_);
        
            vesicle_.Reparam();
        }
    };
};

#endif //_TIMESTEPPER_H_
