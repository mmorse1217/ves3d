#include "ScalarsTest.h"
#include "DeviceCPU.h"

int main(int argc, char* argv[])
{
    DeviceCPU<float> cpu;
    Scalars<float> scfld(cpu);
    ScalarsTest<float> sctest(scfld);
   
    sctest.performAll();
    
    return 0;
}



    


