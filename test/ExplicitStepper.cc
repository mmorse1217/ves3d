#include "DeviceCPU.h"
#include "DataIO.h"
#include "TimeStepper.h"
#include "VesUtil.h"
#include <cmath>
using namespace std;
typedef float T;

int main(int argc, char **argv)
{
//     int n_surfs = 1;
//     int stride = 100;
//     int np = stride * n_surfs;

//     float *x_in = new float[3*np];
//     float *density = new float[3*np];
//     float *vel_out = new float[3*np];

//     for(int ii=0;ii<n_surfs;++ii)
//         for(int jj=0;jj<stride;++jj)
//         {
//             int idx = 3*ii*stride + jj;
//             int m = ii*stride + jj;
//             x_in[idx           ] = cos((2*M_PI*m)/np);
//             x_in[idx +   stride] = sin((2*M_PI*m)/np);
//             x_in[idx + 2*stride] = 0;
            
//             density[idx           ] = cos((2*M_PI*m)/np);
//             density[idx +   stride] = sin((2*M_PI*m)/np);
//             density[idx + 2*stride] = 0;

//         }
    
//     DirectInteraction(x_in, x_in, stride, n_surfs, vel_out);

//     for(int ii=0;ii<n_surfs;++ii)
//         for(int jj=0;jj<stride;++jj)
//         {
//             int idx = 3*ii*stride + jj;
//             //printf("%2.4f\t %2.4f\t %2.4f\t\n", vel_out[idx],vel_out[idx + stride],vel_out[idx + 2*stride]);
//             cout<<vel_out[idx] * vel_out[idx] + vel_out[idx + stride] * vel_out[idx + stride]<<endl;
//         }
    
//     delete[] x_in;
//     delete[] density;
//     delete[] vel_out;

//     return 0;

    //    int n_lattice = 4;

    //Surface parameters
    SurfaceParams<T> par;
    par.p_ = 12;
    par.n_surfs_ = 2;//n_lattice * n_lattice * n_lattice;
    par.kappa_ = 1e-2;
    par.filter_freq_ = 8;
    par.rep_ts_ = 5e-2;
    par.rep_max_vel_ = 5e-4;
    par.rep_iter_max_ = 100;
    par.rep_up_freq_ = 24;
    par.rep_filter_freq_ = 4;

    //Time stepping parameters
    T ts(1e-2); //1e-2 works for single, 
    int n_steps(20000);
    
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
    DataIO<T> myIO(cpu_device,"Positions7.txt", data_length);
    char fname[100];
    sprintf(fname,"../precomputed/biconcave_ra95_%u",par.p_);
    //sprintf(fname,"../precomputed/dumbbell_cart%u_single.txt",par.p_);
    myIO.ReadData(fname, one_vec_length, vesicle.x_.data_);

    //populate copies
    T centers[6] ={0, .5, 0,
                   3,0, 0};
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
