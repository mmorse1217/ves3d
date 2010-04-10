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
    
    //Length and other derived parameters
    int one_vec_length = 6 * surf_par.p_ * (surf_par.p_+1);
    int data_length = one_vec_length * surf_par.n_surfs_;

    //Time stepping parameters
    T ts(time_par.ts_); //1e-2 works for single, 
    int n_steps(time_par.n_steps_);
    
    //Background flow
    ShearFlow<T> flow_field;
    flow_field.shear_rate_ = time_par.shear_rate_;

    //Setting up the device
    DeviceCPU<T> cpu_device;

    //IO 
    DataIO<T> myIO(cpu_device,file_list.simulation_out_file_, data_length);

    //Reading matrices
    bool readFromFile = true;
    OperatorsMats<T> mats(myIO,surf_par.p_, surf_par.rep_up_freq_, readFromFile);

    //Initializing the device
    cpu_device.InitializeSHT(mats);

    //Setting up the surface
    Surface<T> vesicle(cpu_device, surf_par,mats);
    
    //Reading initial shape to the first vesicle
    char fname[200]; 
    file_list.initial_shape_file_ = file_list.initial_shape_file_ + "%u";
    sprintf(fname,file_list.initial_shape_file_.c_str(),surf_par.p_);
    myIO.ReadData(fname, one_vec_length, vesicle.x_.data_);
    
    for(int ii=1;ii<surf_par.n_surfs_;++ii)
      memcpy(vesicle.x_.data_ + ii*one_vec_length, vesicle.x_.data_, one_vec_length *sizeof(T));

    vesicle.UpdateAll();

    //Time stepper
    TimeStepper<T> exp_stepper(ts, n_steps, vesicle, myIO, flow_field, 
			       mats.quad_weights_p_up_, NULL);//&DirectInteraction);

    exp_stepper.saveData = false;
    exp_stepper.verbose = true;
    exp_stepper.userMonitor =NULL;
    exp_stepper.user = (void*) vesicle.work_arr;


#ifdef PROFILING_LITE
    double ss = get_seconds();
#endif
    exp_stepper.EvolveInTime();

#ifdef PROFILING_LITE
    ss = get_seconds()-ss ;
    cout<<" . EvolveInTime (sec) : "<<ss<<endl;
#endif

    
    cout<<" . Total Flops : "<<Logger::GetGFlops()<< "GFlops."<<endl;
}
