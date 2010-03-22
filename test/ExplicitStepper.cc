#include "DeviceCPU.h"
#include "DataIO.h"
#include "TimeStepper.h"
#include "VesUtil.h"

using namespace std;
typedef float T;

int main(int argc, char **argv)
{
    //Surface parameters
    SurfaceParams<T> par;
    par.p_ = 12;
    par.n_surfs_ = 100;
    par.kappa_ = 1e-2;
    par.filter_freq_ = 8;
    par.rep_ts_ = 0.1;
    par.rep_max_vel_ = 1e-2;
    par.rep_iter_max_ = 10;
    par.rep_up_freq_ = 24;
    par.rep_filter_freq_ = 4;

    //Time stepping parameters
    T ts(1);
    int n_steps(5);
    
    //Background flow
    T shear_rate = .1;
    ShearFlow<T> flow_field;
    flow_field.shear_rate_ = shear_rate;
    
    //Setting up the device
    DeviceCPU<T> cpu_device;

    //Initializing the device
    cpu_device.InitializeSHT(par.p_, par.rep_up_freq_);

    //Length and other derived parameters
    int one_vec_length = 6 * par.p_ * (par.p_+1);
    int data_length = one_vec_length * par.n_surfs_;
 
    //Setting up the surface
    Surface<T> vesicle(cpu_device, par);
    
    //Reading initial positions
    DataIO<T> myIO;
    char fname[100];
    sprintf(fname,"../data/dumbbell_cart%u_single.txt",par.p_);
    myIO.ReadData(fname, one_vec_length, vesicle.x_.data_);

    //populate copies
    for(int ii=1;ii<par.n_surfs_;ii++)
        for(int idx=0;idx<one_vec_length;idx++)
            vesicle.x_.data_[ii*one_vec_length + idx] = 3*ii + vesicle.x_.data_[idx];
    vesicle.UpdateAll();

    //Time-stepper
#ifdef PROFILING
    double ss = get_seconds();
#endif
    TimeStepper<T> exp_stepper(ts, n_steps, vesicle, flow_field, &DirectInteraction);
    exp_stepper.EvolveInTime();

#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<" The whole simulation (sec) : "<<ss<<endl;
#endif

    // save the final result
    sprintf(fname,"../data/cart%u_final.txt",par.p_);
    myIO.WriteData(fname,vesicle.x_.GetDataLength(),vesicle.x_.data_);
}
