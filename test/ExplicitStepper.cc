#include "DeviceCPU.h"
#include "DataIO.h"
#include "TimeStepper.h"
#include "VesUtil.h"

using namespace std;
typedef float T;

int main(int argc, char **argv)
{
    //    int n_lattice = 4;

    //Surface parameters
    SurfaceParams<T> par;
    par.p_ = 12;
    par.n_surfs_ = 2;//n_lattice * n_lattice * n_lattice;
    par.kappa_ = 1e-2;
    par.filter_freq_ = 8;
    par.rep_ts_ = 1e-1;
    par.rep_max_vel_ = 5e-4;
    par.rep_iter_max_ = 50;
    par.rep_up_freq_ = 24;
    par.rep_filter_freq_ = 4;

    //Time stepping parameters
    T ts(2e-2); //1e-2 works for single, 
    int n_steps(2000);
    
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
    
    //Reading initial shape to the first vesicle
    DataIO<T> myIO(cpu_device,"Positions.txt", 2*data_length);
    char fname[100];
    sprintf(fname,"../precomputed/biconcave_ra95_%u",par.p_);
    //sprintf(fname,"../precomputed/dumbbell_cart%u_single.txt",par.p_);
    myIO.ReadData(fname, one_vec_length, vesicle.x_.data_);

    //populate copies
    T centers[6] ={0, 0, 0,
                   2.2,0, 0};
//    T *centers = (T*) malloc(3*par.n_surfs_ * sizeof(T));

//     for(int ii=0;ii<n_lattice;ii++)
//         for(int jj=0;jj<n_lattice;jj++)
//             for(int kk=0;kk<n_lattice;kk++)
//             {
//                 int idx = ii * n_lattice * n_lattice + jj * n_lattice + kk;
//                 centers[3*idx  ] = 1.2 * kk;
//                 centers[3*idx+1] = 2.2 * jj;
//                 centers[3*idx+2] = 2.2 * ii - ((n_lattice-1) * 2.2/2);
//             }
    
    vesicle.Populate(centers);
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

}
