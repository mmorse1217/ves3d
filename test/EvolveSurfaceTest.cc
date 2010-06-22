#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
#include "EvolveSurface.h"
#include "OperatorsMats.h"

extern const Device<CPU> the_device(0);

typedef float real;
typedef float fmm_value_type;

int main(int argc, char **argv)
{
    typedef containers::Scalars<real, CPU, the_device> Sca;
    typedef containers::Vectors<real, CPU, the_device> Vec;
    typedef Surface<Sca,Vec> Sur;
    typedef Parameters<real> Par;
    typedef VesInteraction<fmm_value_type> Interaction;
    
    // Setting the parameters
    Par::getInstanceModifiable().n_surfs = 1;   
    Par::getInstanceModifiable().ts = .1;    
    Par::getInstanceModifiable().time_horizon = 1;
    Par::getInstanceModifiable().inner_solver_maxit = 15;    
    Par::getInstanceModifiable().bg_flow_param = 0.1;    
    cout<<Par::getInstance()<<endl;

    //IO
    DataIO<real, CPU> myIO(the_device);
    
    //Initializing vesicle positions from text file
    Vec x0(Par::getInstance().n_surfs,
        Par::getInstance().sh_order);
    
    //reading the prototype form file
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0.getSubN(1)-x0.begin(), x0.begin());
    
    //Making centers and populating the prototype
    Vec cntrs(x0.getNumSubs(), 0, make_pair(1,1));
    real cntrs_host[] = {-5, 0,  1,
                         5, 0, -1};
   
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    
    Populate(x0, cntrs);

    // The interaction class
    Interaction interaction(NULL);//&StokesAlltoAll);

    //Reading operators from file
    bool readFromFile = true;
    OperatorsMats<real> mats(myIO, readFromFile);

    //Making the surface, and time stepper
    Sur S(x0, mats);
    EvolveSurface<Sur, Interaction> Es(mats);
    Es(S, interaction);

    myIO.WriteData("EvolveSurf.out", S.getPosition().size(), 
        S.getPosition().begin());

    PROFILEREPORT(SortTime);
}
