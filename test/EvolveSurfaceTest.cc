#include "EvolveSurface.h"
#include "Device.h"
#include "Scalars.h"
#include "Vectors.h"
#include "Surface.h"
#include "Parameters.h"

extern const Device<CPU> the_cpu_dev(0);
extern const Device<GPU> the_gpu_dev(0);

typedef double real;

int main(int argc, char **argv)
{
    typedef containers::Scalars<real, CPU, the_cpu_dev> Sca;
    typedef containers::Vectors<real, CPU, the_cpu_dev> Vec;
    typedef Surface<Sca,Vec> Sur;
    typedef Parameters<real> Par;
    
    // Setting the parameters
    Par::getInstanceModifiable().n_surfs = 2;   
    Par::getInstanceModifiable().n_steps = 40;
    Par::getInstanceModifiable().time_horizon = 20;
    Par::getInstanceModifiable().ts = .5;    
    
    cout<<Par::getInstance()<<endl;
    
    //IO
    DataIO<Sca::value_type, CPU> myIO(the_cpu_dev);
    
    //Initializing vesicle positions from text file
    Vec x0(Par::getInstance().n_surfs,
        Par::getInstance().sh_order);
    
    //reading the prototype form file
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0.getSubN(1)-x0.begin(), x0.begin());
    
    //Making centers and populating the prototype
    Vec cntrs(x0.getNumSubs(), 0, make_pair(1,1));
    real cntrs_host[] = {0, 0,  1,
                         0, 2, -1};
    
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    
    Populate(x0, cntrs);

    //Making the surface, and time stepper
    Sur S(x0);
    EvolveSurface<Sur> Es;
    Es(S);

    myIO.WriteData("EvolveSurf.txt", S.getPosition().size(), 
        S.getPosition().begin());

    PROFILEREPORT(SortTime);
}
