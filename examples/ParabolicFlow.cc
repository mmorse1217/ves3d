#include "EvolveSurface.h"
#include "CPUKernels.h"

#define DT CPU 
typedef double real;
extern const Device<DT> the_device(0);

int main(int argc, char **argv)
{
    COUT("\n\n ================================"<<
           "\n  Single Vesicle, parabolic flow: "<<
           "\n ================================"<<endl);
    
    typedef Parameters<real> Par_t;
    typedef EvolveSurface<real, DT, the_device> Evolve_t;
    typedef Evolve_t::Sca_t Sca_t;
    typedef Evolve_t::Vec_t Vec_t;

    // Setting the parameters
    int n_layer = 1;
    int n0 = 1;

    Par_t sim_par;
    sim_par.n_surfs = n_layer * n0;
    sim_par.ts = .1;    
    sim_par.time_horizon = 20;
    sim_par.bending_modulus = 1e-2;
    sim_par.rep_maxit = 200;
    sim_par.save_data = true;    
    sim_par.save_stride = .2;
    sim_par.save_file_name = "ParabolicFlow2.txt";
    sim_par.rep_tol = 1e-6;

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
    int length = DIM * n0;
    Array<real, DT, the_device> cntrs(length);
    myIO.ReadData("precomputed/pepperoni.txt", cntrs, 0, length);
    axpy(.5, cntrs, cntrs); //scaling
    cntrs.resize(length * n_layer);
    for(int ii=1; ii<n_layer; ++ii)
        cntrs.getDevice().Memcpy(cntrs.begin() + ii*length,
            cntrs.begin(), length * sizeof(real), MemcpyDeviceToDevice);
    
    Array<real, DT, the_device> perturb(cntrs.size());
    perturb.fillRand();
    axpy(1.0,cntrs, perturb, cntrs);
    //cntrs.begin()[2] = 1;
    Populate(x0, cntrs);

    //Reading Operators From File
    bool readFromFile = true;
    Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    real radius(3), center_vel(1);
    ParabolicFlow<Sca_t, Vec_t> vInf(radius, center_vel);
                
    Evolve_t Es(sim_par, Mats, x0, &vInf, &StokesAlltoAll);
        
    CLEARERRORHIST();
    PROFILESTART(); 
    QC( Es.Evolve() );    
    PROFILEEND("",0);
    
    PRINTERRORLOG();
    PROFILEREPORT(SortFlop);
}
