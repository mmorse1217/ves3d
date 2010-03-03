#include "VectorsTest.h"
#include "DeviceCPU.h"

int main(int argc, char* argv[])
{
    DeviceCPU<float> cpu;
    VectorsTest<float> vectest(cpu);
    vectest.performAll();
    
    return 0;
}



    


