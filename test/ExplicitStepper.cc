#include "TimeStepper.h"
#include "DataIO.h"
#include "DeviceCPU.h"

using namespace std;
typedef float T;

int main(int argc, char **argv)
{

    int p(12), num_vesicles(1);
    T ts(1e-3);
    int n_steps(10);

    //Length and other derived parameters
    int data_length = 6 * p * (p+1) * num_vesicles;

    //Setting up the device
    DeviceCPU<T> cpu_device;
    DataIO<T> myIO;

    cpu_device.InitializeSHT(p, "../data/legTrans12_single.txt",
        "../data/legTransInv12_single.txt",
        "../data/d1legTrans12_single.txt",
        "../data/d2legTrans12_single.txt");
    
    //Setting up the surface
    Surface<T> vesicle(cpu_device, p, num_vesicles);
    
    //Reading initial positions
    myIO.ReadData("../data/dumbbell_cart12_single.txt",data_length,
        vesicle.x_.data_);
    vesicle.UpdateProps();
    //populate copies
    
    //Time-stepper
    TimeStepper<T> exp_stepper(ts, n_steps, vesicle);
    exp_stepper.EvolveInTime();

    // save the final result
    myIO.WriteData("../data/cart12_final.txt",data_length,vesicle.x_.data_);
}
