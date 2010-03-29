#include "DeviceCPU.h"
#include "DataIO.h"
#include "TimeStepper.h"
#include "VesUtil.h"
#include "Logger.h"
#include "Monitor.h"
#include <cmath>

using namespace std;
typedef float T;

int main(int argc, char *argv[])
{
    //Reading parameters from file
    AllParams<T> all_par;   
    all_par = ParserLite<AllParams<T> >("ExplicitStepper.conf");

    SurfaceParams<T> surf_par = all_par.s_par;
    TimeStepperParams<T> time_par = all_par.t_par;
    FileList file_list = all_par.f_list;
    
    //Time stepping parameters
    T ts(time_par.ts_); //1e-2 works for single, 
    int n_steps(time_par.n_steps_);
    
    //Background flow
    ShearFlow<T> flow_field;
    flow_field.shear_rate_ = time_par.shear_rate_;
    
    //Setting up the device
    DeviceCPU<T> cpu_device;
    bool readFromFile = true;
    OperatorsMats<T> mats(surf_par.p_, surf_par.rep_up_freq_, readFromFile);

    //Initializing the device
    cpu_device.InitializeSHT(mats);

    //Length and other derived parameters
    int one_vec_length = 6 * surf_par.p_ * (surf_par.p_+1);
    int data_length = one_vec_length * surf_par.n_surfs_;
 
    //Setting up the surface
    Surface<T> vesicle(cpu_device, surf_par,mats);
    
    //Reading initial shape to the first vesicle
    DataIO<T> myIO(cpu_device,file_list.simulation_out_file_, data_length);

    char fname[200]; 
    file_list.initial_shape_file_ += "%u";
    sprintf(fname,file_list.initial_shape_file_.c_str(),surf_par.p_);
    myIO.ReadData(fname, data_length, vesicle.x_.data_);

    vesicle.UpdateAll();

    //The monitor
    MntrOpts mntr_opts;
    
    mntr_opts.save_centers_ = false;
    mntr_opts.save_shapes_ = false;
    mntr_opts.save_freq_ = 1;
    mntr_opts.area_inc_fac_ = 2;

    Monitor<T> mntr(mntr_opts, myIO);
    
    //Time stepper
    TimeStepper<T> exp_stepper(ts, n_steps, vesicle, myIO, flow_field, 
        mats.quad_weights_p_up_, mntr, &DirectInteraction);

    exp_stepper.EvolveInTime();
    
    cout<<"Total Flops : "<<Logger::GetGFlops()<< "GFlops."<<endl;
}
