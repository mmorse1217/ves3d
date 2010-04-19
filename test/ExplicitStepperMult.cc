#include "DeviceGPU.h"
#include "DeviceCPU.h"
#include "DataIO.h"
#include "TimeStepperMult.h"
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
    all_par = ParserLite<AllParams<T> >("ExplicitStepperAll.conf");
    
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

    //Vesicle vector
    Surface<T> **vesicles = new Surface<T>*[2];

    DeviceCPU<T> cpu_device;

    //IO classes
    DataIO<T> cpuIO(cpu_device,file_list.simulation_out_file_ + "cpu", 5*data_length);

    //Reading operators from file
    bool readFromFile = true;
    OperatorsMats<T> cpu_mats(cpuIO, surf_par.p_, surf_par.rep_up_freq_, readFromFile);

    //Initializing the device
    cpu_device.InitializeSHT(cpu_mats);    
    
    //Background flow
    ShearFlow<T> flow_field;
    flow_field.shear_rate_ = time_par.shear_rate_;

    //Setting up the surface
    vesicles[0] = new Surface<T>(cpu_device, surf_par,cpu_mats);

    //Reading the base shape
    char fname[200]; 
    file_list.initial_shape_file_ += "%u";
    sprintf(fname,file_list.initial_shape_file_.c_str(),surf_par.p_);
    cpuIO.ReadData(fname, one_ves_length, vesicles[0]->x_.data_);
    
    //Making centers
    T *cnts_host = new T[3 * surf_par.n_surfs_];
    T h[3] = {5,5,5};
    T zero[3] = {0, 0, 0};
    for(int k=0; k<nd; k++){
        for(int j=0; j<nd;j++){ 
            for(int i=0; i<nd; i++){
                unsigned int idx = 3*k*nd*nd + 3*j*nd + 3*i;
                cnts_host[idx   ] = zero[0] + i*h[0];
                cnts_host[idx +1] = zero[1] + j*h[1];
                cnts_host[idx +2] = zero[2] + k*h[2];
            }
        }
    }
    
    //Populate cpu
    vesicles[0]->Populate(cnts_host);
    delete[] cnts_host;
    
    //Updating
    vesicles[0]->UpdateAll();
    cout<<" - Populated and updated on cpu."<<endl;
            
    int device_id= 0;
    DeviceGPU<T> gpu_device(device_id);
    
    string fstring = file_list.simulation_out_file_ + "gpu_";
    DataIO<T> gpuIO(gpu_device, fstring.c_str(), 2*data_length);
    
    OperatorsMats<T> gpu_mats(gpuIO, surf_par.p_, surf_par.rep_up_freq_, false);
    gpu_device.Memcpy(gpu_mats.data_, cpu_mats.data_, cpu_mats.GetDataLength(),MemcpyHostToDevice);
    cout<<" - MATS loaded and copied to gpu "<<endl;
    
    gpu_device.InitializeSHT(gpu_mats);
    cout<<" - SHT intialized on gpu "<<endl;
    
    vesicles[1] = new Surface<T>(gpu_device, surf_par,gpu_mats);
    cout<<" - Surface built on gpu "<<endl;
    
    gpu_device.Memcpy(vesicles[1]->x_.data_, vesicles[0]->x_.data_, one_ves_length, MemcpyHostToDevice);
    
    //Making centers
    cnts_host = new T[3 * surf_par.n_surfs_];
    zero[0] = nd * h[0];
    for(int k=0; k<nd; k++){
        for(int j=0; j<nd;j++){ 
            for(int i=0; i<nd; i++){
                unsigned int idx = 3*k*nd*nd + 3*j*nd + 3*i;
                cnts_host[idx   ] = zero[0] + i*h[0];
                cnts_host[idx +1] = zero[1] + j*h[1];
                cnts_host[idx +2] = zero[2] + k*h[2];
            }
        }
    }
    delete[] cnts_host;

    //Populate gpu
    vesicles[1]->Populate(cnts_host); //NOTE THAT CENTERS ARE ALWAYS ON THE HOST
    
    vesicles[1]->UpdateAll();
    cout<<" - Populated and updated."<<endl;

    int num_ves = 2;
    {//The section is required to assure that the TimeStepperMult is
     //destroyed before the vesicles
        TimeStepperMult<T> stepper(ts, n_steps, vesicles, num_ves,  cpuIO, flow_field, 
            cpu_mats.quad_weights_p_up_, NULL);
        
        stepper.saveData = false;
        stepper.verbose = true;
        stepper.userMonitor =NULL;
        stepper.EvolveInTime();
    }

    delete vesicles[0];
    delete vesicles[1];
    delete[] vesicles;

    return 0;
}
