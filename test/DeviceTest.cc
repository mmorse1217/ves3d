//#include "DeviceGPU.cu"
#include "DeviceCPU.h"
#include "DeviceTest.h"

int main(int argc, char* argv[])
{
    DeviceCPU<float> cpu_f;
    cpu_f.InitializeSHT(12,24);
    DeviceTest<float> dvt_f(cpu_f);
    dvt_f.PerformAll();

    return 0;
}
