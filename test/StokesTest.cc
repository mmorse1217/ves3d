#include "Vectors.h"

//#include "InterfacialVelocity.h"

typedef float real;

extern const Device<CPU> the_cpu_device(0);

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Stokes Test:"
        <<"\n ==============================\n");
    sleep(1);

    typedef Scalars<real, CPU, the_cpu_device> ScaCPU_t;
    typedef Vectors<real, CPU, the_cpu_device> VecCPU_t;
    
    return 0;
}
