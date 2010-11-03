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
    sim_par.ts = 14;    
    sim_par.time_horizon = 2800;
    sim_par.rep_maxit = 20;
    sim_par.save_data = false;    
    sim_par.save_stride = 0;
    sim_par.save_file_name = "EvolveSurf.out";
    
    sim_par.scheme = SemiImplicit;
    sim_par.bg_flow_param = 0;

    int p = 48;
    sim_par.sh_order = p;
    sim_par.filter_freq = 2*p/3;
    sim_par.rep_up_freq = 2*p;
    sim_par.rep_filter_freq = p/3;
    sim_par.singular_stokes = ViaSpHarm;

    sim_par.position_solver_iter = 2*p;
    sim_par.tension_solver_iter = 2*p;

    COUT(sim_par<<endl);
    remove(sim_par.save_file_name.c_str());

    //Initial vesicle positions 
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file  
    DataIO myIO(sim_par.save_file_name);
    char fname[300];
    string prec = (typeid(real) == typeid(float)) ? "float" : "double"; 
    sprintf(fname,"precomputed/dumbbell_%u_%s.txt",sim_par.sh_order,prec.c_str());
    myIO.ReadData(fname, x0, 0, x0.getSubLength());
    
    //Reading Operators From File
    bool readFromFile = true;
    Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);
                
    Evolve_t Es(sim_par, Mats, x0, &vInf, NULL);
        
    CLEARERRORHIST();
    PROFILESTART(); 
    QC( Es.Evolve() );    
    PROFILEEND("",0);
    
    PRINTERRORLOG();
    PROFILEREPORT(SortFlop);
    myIO.Append(Es.S_->getPosition());
}
