#include "DeviceGPU.h"
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
    all_par = ParserLite<AllParams<T> >("ExplicitStepperGPU.conf");

    SurfaceParams<T> surf_par = all_par.s_par;
    TimeStepperParams<T> time_par = all_par.t_par;
    FileList file_list = all_par.f_list;

    //Length and other derived parameters
    int nd = surf_par.n_surfs_; //number of vesicles per dimension
    surf_par.n_surfs_ = surf_par.n_surfs_ * surf_par.n_surfs_ * surf_par.n_surfs_;
    int one_ves_length = 6 * surf_par.p_ * (surf_par.p_+1);
    int data_length    = one_ves_length * surf_par.n_surfs_;
    
    //Time stepping parameters
    T ts(time_par.ts_);
    int n_steps(time_par.n_steps_);
    
    //Setting up the device
    int device_id=0;
    DeviceGPU<T> gpu_device(device_id);
    
    //IO classes
    DataIO<T> gpuIO(gpu_device,file_list.simulation_out_file_ + "gpu", 10*data_length);

    //Reading operators from file
    bool readFromFile = true;
    OperatorsMats<T> mats_gpu(gpuIO, surf_par.p_, surf_par.rep_up_freq_, readFromFile);
    cout<<" - MATS loaded"<<endl;

    //Initializing the device
    gpu_device.InitializeSHT(mats_gpu);
    cout<<" - SHT intialized"<<endl;

    //Background flow
    ParabolicFlow<T> flow_field;
    flow_field.R = 10;
    flow_field.U = 3;

    //Setting up the surface
    Surface<T> vesicle_gpu(gpu_device, surf_par,mats_gpu);
    cout<<" - Surface built"<<endl;

    char fname[200]; 
    file_list.initial_shape_file_ += "%u";
    sprintf(fname,file_list.initial_shape_file_.c_str(),surf_par.p_);
    gpuIO.ReadData(fname, one_ves_length, vesicle_gpu.x_.data_);

    //Making centers
    T *cnts_host = new T[3 * surf_par.n_surfs_];

    //Populate
    vesicle_gpu.Populate(cnts_host); //NOTE THAT CENTERS ARE ALWAYS ON THE HOST
    vesicle_gpu.UpdateAll();
    cout<<" - Populated"<<endl;
    delete[] cnts_host;
    
    //Time stepper -- no interaction
    TimeStepper<T> stepper_gpu(ts, n_steps, vesicle_gpu, gpuIO, flow_field, 
        mats_gpu.quad_weights_p_up_, &DirectInteraction);//NULL);
    
    stepper_gpu.saveData = false;
    stepper_gpu.verbose = true;
    stepper_gpu.userMonitor =NULL;
    stepper_gpu.user = (void*) vesicle_gpu.work_arr;
    
    //Evolve
    stepper_gpu.EvolveInTime();
    
    cout<<"Total Flops : "<<Logger::GetGFlops()<< "GFlops."<<endl;

    gpuIO.Append(vesicle_gpu.x_.data_, vesicle_gpu.x_.GetDataLength());            
    return 0;
}


