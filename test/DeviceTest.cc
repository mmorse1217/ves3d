//#include "DeviceGPU.cu"
#include "DeviceCPU.h"
#include "DeviceTest.h"

int main(int argc, char* argv[])
{
    DeviceCPU<float> cpu_f;
    cpu_f.InitializeSHT(12,"../data/legTrans12", "../data/legTransInv12", "../data/d1legTrans12", "../data/d2legTrans12");
    DeviceTest<float> dvt_f(cpu_f);
    dvt_f.PerformAll();

    float *in;
//     DeviceGPU<float> gpu_f;
//     DeviceTest<float> dvt_f2(gpu_f);
//     dvt_f2.TestMemcpy();

    return 0;
}
