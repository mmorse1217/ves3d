#include "DeviceGPU.h"
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

    DeviceCPU<T> cpu_device;

    //IO classes
    DataIO<T> cpuIO(cpu_device,file_list.simulation_out_file_ + "cpu", 10*data_length);

    //Reading operators from file
    bool readFromFile = true;
    OperatorsMats<T> cpu_mats(cpuIO, surf_par.p_, surf_par.rep_up_freq_, readFromFile);

    //Initializing the device
    cpu_device.InitializeSHT(cpu_mats);    
    
    //Background flow
    ShearFlow<T> flow_field;
    flow_field.shear_rate_ = time_par.shear_rate_;

    //Setting up the surface
    Surface<T> cpu_vesicle(cpu_device, surf_par,cpu_mats);

    //Reading the base shape
    char fname[200]; 
    file_list.initial_shape_file_ += "%u";
    sprintf(fname,file_list.initial_shape_file_.c_str(),surf_par.p_);
    cpuIO.ReadData(fname, one_ves_length, cpu_vesicle.x_.data_);
    
    double ss=get_seconds();
    
    //Setting up the all devices
#pragma omp parallel num_threads(3)
    {
        ///cpu device thread
        if(omp_get_thread_num() == 0)
        {
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
            cpu_vesicle.Populate(cnts_host);
            delete[] cnts_host;

            //Updating
            cpu_vesicle.UpdateAll();
            cout<<" - Populated and updated on cpu."<<endl;
            
            //Time stepper 
            TimeStepper<T> cpu_stepper(ts, n_steps, cpu_vesicle, cpuIO, flow_field, 
                cpu_mats.quad_weights_p_up_, NULL);
            
            cpu_stepper.saveData = false;
            cpu_stepper.verbose = true;
            cpu_stepper.userMonitor =NULL;
            cpu_stepper.user = (void*) cpu_vesicle.work_arr;

            cpu_stepper.EvolveInTime();
            
            cpuIO.Append(cpu_vesicle.x_.data_, cpu_vesicle.x_.GetDataLength());            
        }

        if(omp_get_thread_num() != 0)
        {
            int this_thread = omp_get_thread_num()-1; 
            int device_id=  this_thread;
            DeviceGPU<T> gpu_device(device_id);
                
            char fname[300];
            string fstring = file_list.simulation_out_file_ + "gpu_" + "%u";
            sprintf(fname, fstring.c_str(), this_thread);
            DataIO<T> gpuIO(gpu_device, fname, 2*data_length);
            
            OperatorsMats<T> gpu_mats(gpuIO, surf_par.p_, surf_par.rep_up_freq_, false);
            gpu_device.Memcpy(gpu_mats.data_, cpu_mats.data_, cpu_mats.GetDataLength(),MemcpyHostToDevice);
            cout<<" - MATS loaded and copied to device "<<this_thread<<endl;

            gpu_device.InitializeSHT(gpu_mats);
            cout<<" - SHT intialized on device "<<this_thread<<endl;

            Surface<T> gpu_vesicle(gpu_device, surf_par,gpu_mats);
            cout<<" - Surface built on device "<<this_thread<<endl;
            
            gpu_device.Memcpy(gpu_vesicle.x_.data_, cpu_vesicle.x_.data_, one_ves_length, MemcpyHostToDevice);

            //Making centers
            T *cnts_host = new T[3 * surf_par.n_surfs_];
            T h[3] = {5,5,5};
            T zero[3] = {0, 0, 0};
            zero[0] = omp_get_thread_num() * nd * h[0];
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
            gpu_vesicle.Populate(cnts_host); //NOTE THAT CENTERS ARE ALWAYS ON THE HOST
            
            gpu_vesicle.UpdateAll();
            cout<<" - Populated and updated."<<endl;
            
            TimeStepper<T> gpu_stepper(ts, n_steps, gpu_vesicle, gpuIO, flow_field, 
                gpu_mats.quad_weights_p_up_, NULL);
            
            gpu_stepper.saveData = false;
            gpu_stepper.verbose = true;
            gpu_stepper.userMonitor =NULL;
            gpu_stepper.user = (void*) gpu_vesicle.work_arr;
            
            gpu_stepper.EvolveInTime();
            
            gpuIO.Append(gpu_vesicle.x_.data_, gpu_vesicle.x_.GetDataLength());            
        }
    
    }

ss=get_seconds()-ss;
    
cout<<"Total time :"<<ss<<endl;
cout<<"Total Flops : "<<Logger::GetGFlops()<< "GFlops."<<endl;
return 0;
}
