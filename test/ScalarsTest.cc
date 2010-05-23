#include "ScalarContainerTest.h"
#include "Scalars.h"

extern const Device<CPU> cpu_dev;

int main(int argc, char* argv[])
{
    ScalarsTest< Scalars<float,CPU, cpu_dev> > sctest;
    sctest.PerformAll();

    return 0;
}



    


