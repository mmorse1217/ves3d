#include "ScalarContainerTest.h"
#include "Scalars.h"
#include <unistd.h>  //for sleep()

extern const Device<CPU> cpu_dev(0);
#ifdef GPU_ACTIVE
extern const Device<GPU> gpu_dev(0);
#endif //GPU_ACTIVE

int main(int argc, char* argv[])
{
    COUT("\n ==============================\n"
        <<"  Scalars Test:"
        <<"\n ==============================\n");
    sleep(1);

    ScalarsTest<Scalars<float,CPU, cpu_dev> > sctest_cpu;
    sctest_cpu.PerformAll();

#ifdef GPU_ACTIVE
    ScalarsTest<Scalars<float,GPU, gpu_dev> > sctest_gpu;
    sctest_gpu.PerformAll();
#endif

    return 0;
}
