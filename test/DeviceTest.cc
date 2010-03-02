#include "DeviceCPU.h"
#include "DeviceGPU.cu"
#include "DeviceTest.h"

int main(int argc, char* argv[])
{
//     DeviceTest<float> dvt_f;
//     DeviceCPU<float> cpu_f;
    
//     dvt_f.SetDevice(&cpu_f);
//     dvt_f.PerformAll();

//     DeviceTest<double> dvt_d;
//     DeviceCPU<double> cpu_d;
    
//     dvt_d.SetDevice(&cpu_d);
//     dvt_d.PerformAll();
    
    DeviceTest<float> dvt_f;
    DeviceGPU<float> gpu_f;
    
    dvt_f.SetDevice(&gpu_f);
    dvt_f.PerformAll();

    return 0;
}
