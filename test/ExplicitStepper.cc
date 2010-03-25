#include "DeviceCPU.h"
#include "DataIO.h"
#include "TimeStepper.h"
#include "VesUtil.h"
//#include "Logger.h"
#include <cmath>

using namespace std;
typedef float T;

// unsigned long int Logger::Mflop_count_ = 0;
// const enum LogLevel Logger::the_log_level_ = FlopCount;

int main(int argc, char *argv[])
{
    //Reading parameters from file
    AllParams<T> all_par;   
    all_par = ParserLite<AllParams<T> >("ExplicitStepper.conf");

    SurfaceParams<T> surf_par = all_par.s_par;
    TimeStepperParams<T> time_par = all_par.t_par;
    FileList file_list = all_par.f_list;
    
    cout<<surf_par<<endl;
    cout<<time_par<<endl;
    cout<<file_list<<endl;

    //Time stepping parameters
    T ts(time_par.ts_); //1e-2 works for single, 
    int n_steps(time_par.n_steps_);
    
    //Background flow
    ShearFlow<T> flow_field;
    flow_field.shear_rate_ = time_par.shear_rate_;
    
    //Setting up the device
    DeviceCPU<T> cpu_device;
    
    //Initializing the device
    cpu_device.InitializeSHT(surf_par.p_, surf_par.rep_up_freq_);

    //Length and other derived parameters
    int one_vec_length = 6 * surf_par.p_ * (surf_par.p_+1);
    int data_length = one_vec_length * surf_par.n_surfs_;
 
    //Setting up the surface
    Surface<T> vesicle(cpu_device, surf_par);
    
    //Reading initial shape to the first vesicle
    DataIO<T> myIO(cpu_device,file_list.simulation_out_file_, data_length);

    char fname[200]; 
    file_list.initial_shape_file_ += "%u";
    sprintf(fname,file_list.initial_shape_file_.c_str(),surf_par.p_);
    myIO.ReadData(fname, data_length, vesicle.x_.data_);

    //populate copies
    //    T centers[6] ={0, .5, 0,
    //              3,0, 0};
//    T *centers = (T*) malloc(3*surf_par.n_surfs_ * sizeof(T));

// //     for(int ii=0;ii<n_lattice;ii++)
// //         for(int jj=0;jj<n_lattice;jj++)
// //             for(int kk=0;kk<n_lattice;kk++)
// //             {
// //                 int idx = ii * n_lattice * n_lattice + jj * n_lattice + kk;
// //                 centers[3*idx  ] = 1.2 * kk;
// //                 centers[3*idx+1] = 2.2 * jj;
// //                 centers[3*idx+2] = 2.2 * ii - ((n_lattice-1) * 2.2/2);
// //             }
    
    //vesicle.Populate(centers);
    //free(centers);
    vesicle.UpdateAll();
    
    //Time-stepper
#ifdef PROFILING
    double ss = get_seconds();
#endif

    TimeStepper<T> exp_stepper(ts, n_steps, vesicle, myIO, flow_field, &DirectInteraction);
    exp_stepper.EvolveInTime();
    
#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<" The whole simulation (sec) : "<<ss<<endl;
#endif

    Logger::TearDown();
}
