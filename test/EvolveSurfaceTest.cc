#include "EvolveSurface.h"
#include "Device.h"
#include "Scalars.h"
#include "Vectors.h"
#include "Surface.h"

extern const Device<CPU> the_cpu_dev(0);

typedef double real;

int main(int argc, char **argv)
{
    typedef containers::Scalars<real, CPU, the_cpu_dev> Sca;
    typedef containers::Vectors<real, CPU, the_cpu_dev> Vec;
    typedef Surface<Sca,Vec> Sur;

    int p(12), nSur(1);

    //IO
    DataIO<Sca::value_type, CPU> myIO(the_cpu_dev);
    
    //Initializing vesicle positions from text file
    Vec x0(nSur, p);
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0.size(), x0.begin());
    
    //Making the surface, and time stepper
    Sur S(x0);
    real dt = 1;
    real T = 100 * dt;    
    EvolveSurface<Sur> Es;
    Es(S, T, dt);

    myIO.WriteData("EvolveSurf.txt", S.getPosition().size(), 
        S.getPosition().begin());
}
