#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

#include "Surface.h"
//#include "Vectors.h"

#include <iostream>

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
        for(int idx=0;idx<n_steps_;idx++)
        {
            vesicle_.StokesMatVec(vesicle_.bending_force, velocity_);
            //vesicle_.StokesMatVec(vesicle_.tensile_force_, vel_tension_);
            //Calculate tension
            axpy((T) 0.0, velocity_, vel_tension_, velocity_); //check whether the same input and output works
            axpy(ts_, velocity_, vesicle_.x_, vesicle_.x_);
            vesicle_.UpdateProps();
#ifdef NDEBUG
            std::cout<<idx<<std::endl;
#endif
        }
    };
};

#endif
