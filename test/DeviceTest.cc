#include "DeviceGPU.cu"
#include "DeviceCPU.h"
#include "DeviceTest.h"

int main(int argc, char* argv[])
{
//     DeviceCPU<float> cpu_f;
//     DeviceTest<float> dvt_f(cpu_f);
//     dvt_f.PerformAll();
    
//     DeviceCPU<double> cpu_d;
//     DeviceTest<double> dvt_d(cpu_d);
//     dvt_d.PerformAll();

    DeviceGPU<float> gpu_f;
    DeviceTest<float> dvt_f2(gpu_f);
    dvt_f2.TestMemcpy();

    return 0;
}
