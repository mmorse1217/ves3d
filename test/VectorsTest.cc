#include "VectorsTest.h"
#include "DeviceCPU.h"

int main(int argc, char* argv[])
{
    DeviceCPU<float> cpu;
    Scalars<float> vecfld(cpu);
    VectorsTest<float> vectest(scfld);
   
    vectest.performAll();
    
    return 0;
}



    


