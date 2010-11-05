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
    typedef Evolve_t::Sca_t Sca_t;
    typedef Evolve_t::Vec_t Vec_t;

    // Setting the parameters
    Par_t sim_par;
    sim_par.n_surfs = 2;
    sim_par.ts = .05;    
    sim_par.time_horizon = 100;
    sim_par.bending_modulus = 1e-2;
    sim_par.rep_maxit = 20;
    sim_par.save_data = true;    
    sim_par.save_stride = .5;
    sim_par.save_file_name = "ShearFlow.txt";
    
    sim_par.scheme = SemiImplicit;
    sim_par.singular_stokes = Direct;

    int p = 12;
    sim_par.sh_order = p;
    sim_par.filter_freq = 2*p/3;
    sim_par.rep_up_freq = 2*p;
    sim_par.rep_filter_freq = p/3;
    
    sim_par.position_solver_iter = 2*p;
    sim_par.tension_solver_iter = 2*p;

    COUT(sim_par<<endl);
    remove(sim_par.save_file_name.c_str());

    //Initial vesicle positions 
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file  
    DataIO myIO(sim_par.save_file_name);
    char fname[300];
    sprintf(fname,"precomputed/biconcave_ra65_%u",sim_par.sh_order);
    myIO.ReadData(fname, x0, 0, x0.getSubLength());

    //Reading the Centers and populating
    real cntrs_host[] = {1.2, 0, -.4, -1.2, 0, .4};
    Array<real, DT, the_device> cntrs(DIM * sim_par.n_surfs);
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host, cntrs.size()
        * sizeof(real), MemcpyHostToDevice);
    Populate(x0, cntrs);

    //Reading Operators From File
    bool readFromFile = true;
    Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    real shear_rate(.1);
    ShearFlow<Vec_t> vInf(shear_rate);
                
    Evolve_t Es(sim_par, Mats, x0, &vInf, &StokesAlltoAll);
        
    CLEARERRORHIST();
    PROFILESTART(); 
    QC( Es.Evolve() );    
    PROFILEEND("",0);
    
    PRINTERRORLOG();
    PROFILEREPORT(SortFlop);
}
