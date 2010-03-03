#include "DeviceCPU.h"
#include "DeviceGPU.cu"
#include "DeviceTest.h"
#include "ScalarsTest.h"

int main(int argc, char* argv[])
{
    DeviceCPU<float> cpu_f;
    DeviceTest<float> dvt_f(cpu_f);
    dvt_f.PerformAll();

    DeviceCPU<double> cpu_d;
    DeviceTest<double> dvt_d(cpu_d);
    dvt_d.PerformAll();
    
    ScalarsTest<float> sct;
    sct.performAll();

//     DeviceGPU<float> gpu_f;
//     DeviceTest<float> dvt_f(gpu_f);

//     //dvt_f.SetDevice(&gpu_f);
//     dvt_f.PerformAll();

    return 0;
}
