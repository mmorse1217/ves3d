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
        velocity_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        vel_tension_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_)
    {};

    ~TimeStepper(){};

    void EvolveInTime()
    {
        cout<<vesicle_.params_;
        
        cout<<" ------------------------------------"<<endl;
        cout<<"  Time stepper parameters"<<endl;
        cout<<" ------------------------------------"<<endl;
        cout<<"  ts                : "<<ts_<<endl;
        cout<<"  n_steps           : "<<n_steps_<<endl;
        cout<<" ------------------------------------"<<endl<<endl;

        for(int idx=0;idx<n_steps_;idx++)
        {

            cout<<idx<<endl;
            vesicle_.UpdateAll();
            
            //Calculate stokes
            vesicle_.StokesMatVec(vesicle_.bending_force_, velocity_);
            vesicle_.StokesMatVec(vesicle_.tensile_force_, vel_tension_);
            //Calculate tension
            vesicle_.GetTension(velocity_, vel_tension_, vesicle_.tension_);
            vesicle_.device_.avpw(vesicle_.tension_, vel_tension_.data_,
                velocity_.data_, velocity_.GetFunLength(), vesicle_.params_.n_surfs_, velocity_.data_);

            //Interaction
//             vesicle_.device_.avpw(vesicle_.tension_, vesicle_.tensile_force_.data_,
//                 vesicle_.bending_force_.data_, vesicle_.x_.GetFunLength(), vesicle_.params_.n_surfs_,
//                 vesicle_.tensile_force_.data_);

//             axpb( 1/8/pi, vesicle_.tensile_force_, (T) 0.0, vesicle_.tensile_force_);


//             vesicle_.device_.DirectStokes(vesicle.x_.GetFunLength(), vesicle_.params_.n_surfs_, 
//                 0, vesicle.x_.GetFunLength(), vesicle_.quad_weights_, vesicle_.x_.data_,
//                 vesicle.x_.data_, const T *den, vel_tension_.data_);




            //Advance in time
            axpy(ts_, velocity_, vesicle_.x_, vesicle_.x_);
        
            vesicle_.Reparam();
        }
    };
};

#endif //_TIMESTEPPER_H_
