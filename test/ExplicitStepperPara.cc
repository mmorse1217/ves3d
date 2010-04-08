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
    all_par = ParserLite<AllParams<T> >("ExplicitStepperPara.conf");

    SurfaceParams<T> surf_par = all_par.s_par;
    TimeStepperParams<T> time_par = all_par.t_par;
    FileList file_list = all_par.f_list;

    //Length and other derived parameters
    int one_ves_length = 6 * surf_par.p_ * (surf_par.p_+1);
    int data_length    = one_ves_length * surf_par.n_surfs_;
    
    //Time stepping parameters
    T ts(time_par.ts_);
    int n_steps(time_par.n_steps_);
    
    //Setting up the device
    int device_id=0;
    DeviceGPU<T> device(device_id);
    
    //IO classes
    DataIO<T> IO(device,file_list.simulation_out_file_, data_length);

    //Reading operators from file
    bool readFromFile = true;
    OperatorsMats<T> mats(IO, surf_par.p_, surf_par.rep_up_freq_, readFromFile);
    cout<<" - MATS loaded"<<endl;

    //Initializing the device
    device.InitializeSHT(mats);
    cout<<" - SHT intialized"<<endl;

    //Background flow
    ParabolicFlow<T> flow_field;
    flow_field.R = 30;
    flow_field.U = 5;

    //Setting up the surface
    Surface<T> vesicle(device, surf_par,mats);
    cout<<" - Surface built"<<endl;

    char fname[200]; 
    file_list.initial_shape_file_ += "%u";
    sprintf(fname,file_list.initial_shape_file_.c_str(),surf_par.p_);
    IO.ReadData(fname, one_ves_length, vesicle.x_.data_);

    //Reading centers
    int nv = surf_par.n_surfs_;
    T *cnts = device.Malloc(3 * nv);
    T *cnts_host = new T[3 * nv];
    
    IO.ReadData("precomputed/pepperoni.txt", 3 * nv, cnts);
    device.Memcpy(cnts_host, cnts, 3*nv, MemcpyDeviceToHost);
    
    //Populate
    vesicle.Populate(cnts_host);
    vesicle.UpdateAll();
    cout<<" - Populated"<<endl;
    delete[] cnts_host;
    device.Free(cnts);

     //Time stepper
    TimeStepper<T> stepper(ts, n_steps, vesicle, IO, flow_field, 
        mats.quad_weights_p_up_, &DirectInteraction);
    
    stepper.saveData = false;
    stepper.verbose = true;
    stepper.userMonitor =NULL;
    stepper.user = (void*) vesicle.work_arr;
    
    //Evolve
    stepper.EvolveInTime();
    
    cout<<"Total Flops : "<<Logger::GetGFlops()<< "GFlops."<<endl;

    IO.Append(vesicle.x_.data_, vesicle.x_.GetDataLength());            
    return 0;
}


