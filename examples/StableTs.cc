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
    sim_par.bg_flow_param = 0;

    //Initial vesicle positions 
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
    DataIO myIO;
    char fname[300];
    sprintf(fname,"precomputed/biconcave_ra85_%u",sim_par.sh_order);
    myIO.ReadData(fname, x0, 0, x0.getSubLength());
    
    //Reading Operators From File
    bool readFromFile = true;
    Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);
    
    int nstep = 200;
    real ts = 1024;
    
    bool doIterate(true);
    while ( doIterate )
    {
        ts /= 2;
        if ( ts < 1e-2 )
            break;
        
        sim_par.ts = ts;    
        sim_par.time_horizon = nstep * ts;
        COUT(sim_par<<endl);

        Evolve_t Es(sim_par, Mats, x0, &vInf, NULL);
        
        CLEARERRORHIST();
        PROFILESTART(); 
        QC( Es.Evolve() );    
        PROFILEEND("",0);
        
        doIterate = !ERRORSTATUS();
        PRINTERRORLOG();
    }

    if ( !doIterate )
        COUT("  The stable time step for the following parameters is "
            << ts <<endl<<sim_par);
    else
        COUT("  The stable time step has not reached.");

    PROFILEREPORT(SortFlop);
}



