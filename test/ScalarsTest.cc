#include "ScalarsTest.h"
#include "Scalars.h"
#include "Device.h"

typedef Device<CPU> DevCPU;
typedef Device<GPU> DevGPU;

extern const DevCPU cpu_dev(0);
#ifdef GPU_ACTIVE
extern const DevGPU gpu_dev(0);
#endif //GPU_ACTIVE

int main(int argc, char* argv[])
{
    VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
    COUT("==============================\n"
        <<" Scalars Test:"
        <<"\n==============================");

    ScalarsTest<Scalars<float, DevCPU, cpu_dev> > sctest_cpu;
    sctest_cpu.PerformAll();

#ifdef GPU_ACTIVE
    ScalarsTest<Scalars<float, DevGPU, gpu_dev> > sctest_gpu;
    sctest_gpu.PerformAll();
#endif

    VES3D_FINALIZE();
    return 0;
}
