#include "VectorsTest.h"
#include "Vectors.h"
#include "Device.h"

typedef Device<CPU> DevCPU;
typedef Device<GPU> DevGPU;

extern const DevCPU cpu_dev(0);
#ifdef GPU_ACTIVE
extern const DevGPU gpu_dev(0);
#endif //GPU_ACTIVE

int main(int argc, char* argv[])
{
    COUT("\n ==============================\n"
        <<"  Vectors Test:"
        <<"\n ==============================\n");

    VectorsTest<Vectors<float, DevCPU, cpu_dev> > tcpu;
    tcpu.PerformAll();

#ifdef GPU_ACTIVE
    VectorsTest<Vectors<float, DevCPU, cpu_dev> > tgpu;
    tgpu.PerformAll();
#endif

    return 0;
}
