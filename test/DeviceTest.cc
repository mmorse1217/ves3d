#include "DeviceCPU.h"
#include "DeviceTest.h"

#ifdef GPU_ACTIVE
#include "DeviceGPU.h"
#endif

int main(int argc, char* argv[])
{
  
    bool readFromFile = true;
    //cpu
    cout<<"=========================================="<<endl;
    cout<<"  CPU "<<endl;
    cout<<"=========================================="<<endl;

    DeviceCPU<float> cpu_f;
    DataIO<float> cpuIO(cpu_f," ", 0);
    OperatorsMats<float> cmats(cpuIO, 6, 6, readFromFile);
    cpu_f.InitializeSHT(cmats);
    DeviceTest<float> dvt_f(cpu_f);
    dvt_f.PerformAll();

#ifdef GPU_ACTIVE
    //gpu
    cout<<"=========================================="<<endl;
    cout<<" GPU "<<endl;
    cout<<"=========================================="<<endl;
    DeviceGPU<float> gpu_f(1);
    DataIO<float> gpuIO(gpu_f," ", 0);
    OperatorsMats<float> gmats(gpuIO, 6, 6, readFromFile);
    gpu_f.InitializeSHT(gmats);
    DeviceTest<float> dvt_gf(gpu_f);
    dvt_gf.PerformAll();
#endif

    return 0;
}
