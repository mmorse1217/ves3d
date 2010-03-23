#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

#include "Surface.h"
#include "VesUtil.h"
#include "DataIO.h"

#define PI_8I 1.0/8.0/M_PI

template<typename T> 
class TimeStepper
{
  public:
    T ts_;
    int n_steps_;
    DataIO<T> &fileIO;
    Surface<T> &vesicle_;
    Vectors<T> velocity_, vel_bending_, vel_tension_;

    VelField<T> &bg_flow_;
    void(*Interaction_)(T*, T*, int, int, T*);

    T *quad_weights_;

    TimeStepper(T ts_in, int n_steps_in, Surface<T> &ves_in, 
        DataIO<T> &fileIO_in, VelField<T> &bg_flow_in, 
        void(*Interaction_in)(T*, T*, int, int,  T*)) : 
        ts_(ts_in),
        n_steps_(n_steps_in), 
        fileIO(fileIO_in),
        vesicle_(ves_in),
        velocity_   (vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        vel_bending_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        vel_tension_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        bg_flow_(bg_flow_in),
        Interaction_(Interaction_in)
    {
        int np = 2 * vesicle_.params_.rep_up_freq_ * (vesicle_.params_.rep_up_freq_ + 1);
        quad_weights_ = (T*) malloc(np * sizeof(T));
        char fname[300];
        sprintf(fname,"../precomputed/quad_weights_%u_single.txt",vesicle_.params_.rep_up_freq_);
        fileIO.ReadData(fname, np, quad_weights_);
    };

    ~TimeStepper()
    {
        free(quad_weights_);
    };

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
            avpw(vesicle_.tension_, vesicle_.tensile_force_, vesicle_.bending_force_, vel_tension_);
            //vel_tension_ holds the interaction force
            xvpb(vesicle_.w_,vel_tension_, (T) 0.0, vel_tension_);
            axpb((T) (PI_8I), vel_tension_, (T) 0.0, vel_tension_);

            //up-sampling x_ to V10
            vesicle_.device_.ShAna(vesicle_.x_.data_, vesicle_.work_arr, vesicle_.params_.p_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V10.data_);
            vesicle_.device_.Resample(vesicle_.params_.p_, 3*vesicle_.params_.n_surfs_, vesicle_.params_.rep_up_freq_, 
                vesicle_.V10.data_, vesicle_.shc);
            vesicle_.device_.ShSyn(vesicle_.shc, vesicle_.work_arr, vesicle_.params_.rep_up_freq_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V10.data_);

            //up-sampling vel_tension to V11
            vesicle_.device_.ShAna(vel_tension_.data_, vesicle_.work_arr, vesicle_.params_.p_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V11.data_);
            vesicle_.device_.Resample(vesicle_.params_.p_, 3*vesicle_.params_.n_surfs_, vesicle_.params_.rep_up_freq_, 
                vesicle_.V11.data_, vesicle_.shc);
            vesicle_.device_.ShSyn(vesicle_.shc, vesicle_.work_arr, vesicle_.params_.rep_up_freq_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V11.data_);

            //the self-interaction -- V12
            vesicle_.device_.DirectStokes(vesicle_.V10.GetFunLength(), vesicle_.params_.n_surfs_, 
                0, vesicle_.V10.GetFunLength(), quad_weights_, vesicle_.V10.data_,
                vesicle_.V10.data_, vesicle_.V11.data_, vesicle_.V12.data_);

            ///@todo the multiplication by the quadrature weights is slow
            int stride = vesicle_.V11.GetFunLength();
            for(int ii=0;ii<vesicle_.params_.n_surfs_;++ii)
            {
                vesicle_.device_.xvpb(quad_weights_, vesicle_.V11.data_ + 3*ii*stride,
                    (T) 0.0, stride, 1, vesicle_.V11.data_ + 3*ii*stride);
            }

            //Interaction to V13
            Interaction_(vesicle_.V10.data_, vesicle_.V11.data_,  
                vesicle_.V10.GetFunLength(), vesicle_.params_.n_surfs_, vesicle_.V13.data_);
            
            axpy((T) -1.0, vesicle_.V12, vesicle_.V13, vesicle_.V13);
            
            vesicle_.device_.ShAna(vesicle_.V13.data_, vesicle_.work_arr, vesicle_.params_.rep_up_freq_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V11.data_);
            vesicle_.device_.Resample(vesicle_.params_.rep_up_freq_, 3*vesicle_.params_.n_surfs_, 
                vesicle_.params_.p_, vesicle_.V11.data_, vesicle_.shc);
            vesicle_.device_.ShSyn(vesicle_.shc, vesicle_.work_arr, vesicle_.params_.p_, 
                3*vesicle_.params_.n_surfs_, velocity_.data_);

            //Filtering
            vesicle_.device_.Filter(vesicle_.params_.p_, 3*vesicle_.params_.n_surfs_, 
                velocity_.data_, vesicle_.alpha_p, vesicle_.work_arr, vesicle_.shc, velocity_.data_);

            //DotProduct(velocity_, velocity_, vesicle_.S1);
            //cout<<vesicle_.S1.Max()<<endl;
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

            if((40*(idx+1))%n_steps_ == 0)
                 fileIO.Append(vesicle_.x_.data_, vesicle_.x_.GetDataLength());
        }
    };
};

#endif //_TIMESTEPPER_H_

