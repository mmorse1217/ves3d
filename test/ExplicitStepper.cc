#include "DeviceCPU.h"
#include "DataIO.h"
#include "TimeStepper.h"

using namespace std;
typedef float T;

int main(int argc, char **argv)
{

    int p(12), num_vesicles(1000);
    T ts(1e-3);
    int n_steps(1);

    //Setting up the device
    DeviceCPU<T> cpu_device;
    DataIO<T> myIO;

    cpu_device.InitializeSHT(p, "../data/legTrans12_single.txt",
        "../data/legTransInv12_single.txt",
        "../data/d1legTrans12_single.txt",
        "../data/d2legTrans12_single.txt");

    //Length and other derived parameters
    int one_vec_length = 6 * p * (p+1);
    int data_length = one_vec_length * num_vesicles;

    //Setting up the surface
    Surface<T> vesicle(cpu_device, p, num_vesicles);
    
    //Reading initial positions
    myIO.ReadData("../data/dumbbell_cart12_single.txt", one_vec_length, vesicle.x_.data_);

    //populate copies
    for(int ii=1;ii<num_vesicles;ii++)
        for(int idx=0;idx<one_vec_length;idx++)
            vesicle.x_.data_[ii*one_vec_length + idx] = 3*ii + vesicle.x_.data_[idx];
    vesicle.UpdateAll();
    
    //Time-stepper
    TimeStepper<T> exp_stepper(ts, n_steps, vesicle);
    exp_stepper.EvolveInTime();

    // save the final result
    myIO.WriteData("../data/cart12_final.txt",one_vec_length,vesicle.x_.data_);
}
