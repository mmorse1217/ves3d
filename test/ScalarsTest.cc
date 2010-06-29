#include "ScalarContainerTest.h"
#include "Scalars.h"

extern const Device<CPU> cpu_dev(0);
extern const Device<GPU> gpu_dev(0);

int main(int argc, char* argv[])
{
    ScalarsTest< containers::Scalars<float,CPU, cpu_dev> > sctest_cpu;
    sctest_cpu.PerformAll();

    ScalarsTest< containers::Scalars<float,GPU, gpu_dev> > sctest_gpu;
    sctest_gpu.PerformAll();

    return 0;
}



    


