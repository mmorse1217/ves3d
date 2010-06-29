#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
#include "OperatorsMats.h"
#include "EvolveSurface.h"

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
    Par sim_par;
    sim_par.n_surfs = 1;   
    sim_par.ts = .5;    
    sim_par.time_horizon = 50;
    sim_par.rep_maxit = 20;
    sim_par.bg_flow_param = 0.1;    
    sim_par.save_data = true;    
    sim_par.save_stride = 1;
    sim_par.save_file_name = "EvolveSurf.out";
    COUT(sim_par);
    
    //Cleaning the slate
    remove(sim_par.save_file_name.c_str());

    //IO
    DataIO<real, CPU> myIO(the_device);

    //Initializing vesicle positions from text file
    Vec x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0.getTheDim()*x0.getStride(), x0.begin());
    
    //Making centers and populating the prototype
    if ( sim_par.n_surfs > 1 )
    {
        Vec cntrs(x0.getNumSubs(), 0, make_pair(1,1));
        real cntrs_host[] = {-5, 0,  1,
                             5, 0, -1};
        
        cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
            cntrs.size() * sizeof(real), MemcpyHostToDevice);
        
        Populate(x0, cntrs);
    }

    // The interaction class
    Interaction interaction(NULL);//&StokesAlltoAll);

    //Reading operators from file
    bool readFromFile = true;
    OperatorsMats<real, DataIO<real, CPU> > mats(myIO, readFromFile, sim_par);

    //Making the surface, and time stepper
    Sur S(x0, mats);
    EvolveSurface<Sur, Interaction> Es(mats, sim_par);
    Es(S, interaction);
    
    PROFILEREPORT(SortTime);
}
