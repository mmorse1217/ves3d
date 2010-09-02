#include "EvolveSurface.h"
#include "CPUKernels.h"

#define DT CPU 
typedef double real;
extern const Device<DT> the_device(0);

int main(int argc, char **argv)
{
    COUT("\n\n ============================="<<
           "\n  Single Vesicle, shear flow: "<<
           "\n ============================="<<endl);
    
    typedef Parameters<real> Par_t;
    typedef EvolveSurface<real, DT, the_device> Evolve_t;
    typedef Evolve_t::Vec_t Vec_t;

    // Setting the parameters
    Par_t sim_par;
    sim_par.n_surfs = 1;
    sim_par.rep_maxit = 20;
    sim_par.save_data = false;    

    sim_par.scheme = Explicit;
    sim_par.bg_flow_param = 0.1;

    //Initial vesicle positions 
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
    DataIO myIO;
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0, 0, x0.getSubLength());
    
    //Reading Operators From File
    bool readFromFile = true;
    Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);
   
    enum MonitorReturn es_status(AreaErrorLarge);
    
    int nstep = 200;
    real ts = 1000;

    while ( es_status != TimeHorizonReached )
    {
        ts /= 2;
        if ( ts < 1e-2 )
            break;

        sim_par.ts = ts;    
        sim_par.time_horizon = nstep * ts;
        COUT(sim_par);

        Evolve_t Es(sim_par, Mats, x0, &vInf, NULL);
        
        PROFILESTART();    
        es_status = Es.Evolve();    
        PROFILEEND("",0);
    }
    PROFILEREPORT(SortFlop);        

}



