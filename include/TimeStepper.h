#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

#include "Surface.h"
#include "VesUtil.h"

#define PI_8I 1.0/8.0/M_PI

template<typename T> 
class TimeStepper
{
  public:
    T ts_;
    int n_steps_;
    Surface<T> &vesicle_;
    Vectors<T> velocity_, vel_bending_, vel_tension_;
    VelField<T> &bg_flow_;
    void (*Interaction_) (T *x_in, T *density_in, T *vel_out);
    

    TimeStepper(T ts_in, int n_steps_in, Surface<T> &ves_in, 
        VelField<T> &bg_flow_in,
        void (*Interaction_in) (T *x_in, T *density_in, T *vel_out) ) : 
        ts_(ts_in),
        n_steps_(n_steps_in), 
        vesicle_(ves_in),
        velocity_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        vel_bending_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        vel_tension_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        bg_flow_(bg_flow_in),
        Interaction_(Interaction_in)
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

            //Interaction
            vesicle_.device_.avpw(vesicle_.tension_, vesicle_.tensile_force_.data_,
                vesicle_.bending_force_.data_, vesicle_.x_.GetFunLength(), vesicle_.params_.n_surfs_,
                vel_tension_.data_); //vel_tension_ holds the interaction force
            
            axpb((T) (PI_8I), vel_tension_, (T) 0.0, vel_tension_);
            
            Interaction_(vesicle_.x_.data_, vel_tension_.data_, vel_bending_.data_);// vel_bending now holds the interaction velocity
            
            vesicle_.device_.DirectStokes(vesicle_.x_.GetFunLength(), vesicle_.params_.n_surfs_, 
                0, vesicle_.x_.GetFunLength(), vesicle_.quad_weights_, vesicle_.x_.data_,
                vesicle_.x_.data_, vel_tension_.data_, velocity_.data_);
            
             axpy((T) -1.0, vel_tension_, velocity_, velocity_);
            
            //Background flow
            bg_flow_.GetVel(vesicle_.x_, vel_tension_);
            axpy((T) 1.0, vel_tension_, velocity_, velocity_);
            
            //Calculate stokes
            vesicle_.StokesMatVec(vesicle_.bending_force_, vel_bending_);
            axpy((T) 1.0, vel_bending_, velocity_, velocity_);

            vesicle_.StokesMatVec(vesicle_.tensile_force_, vel_tension_);

            //Calculate tension
            vesicle_.GetTension(velocity_, vel_tension_, vesicle_.tension_);
            vesicle_.device_.avpw(vesicle_.tension_, vel_tension_.data_,
                velocity_.data_, velocity_.GetFunLength(), vesicle_.params_.n_surfs_, velocity_.data_);

            //Advance in time
            axpy(ts_, velocity_, vesicle_.x_, vesicle_.x_);
        
            vesicle_.Reparam();
        }
    };
};

#endif //_TIMESTEPPER_H_

