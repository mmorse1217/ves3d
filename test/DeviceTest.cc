#include "DeviceGPU.cu"
#include "DeviceCPU.h"
#include "DeviceTest.h"

int main(int argc, char* argv[])
{
//     DeviceCPU<float> cpu_f;
//     cpu_f.InitializeSHT(6,6);
//     //cpu_f.InitializeSHT(6,12); // before uncommenting make sure you have precomputed files.
//     DeviceTest<float> dvt_f(cpu_f);
//     dvt_f.PerformAll();

    DeviceGPU<float> gpu_f;
    DeviceTest<float> dvt_f(gpu_f);
    dvt_f.PerformAll();
    
    return 0;
}
