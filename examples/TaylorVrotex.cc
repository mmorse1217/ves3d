#include "EvolveSurface.h"
#include "CPUKernels.h"

#define DT CPU 
typedef double real;

extern const Device<DT> the_device(0);

int main(int argc, char **argv)
{
    COUT("\n\n ========================\n  Taylor Vortex: "
        <<"\n ========================"<<endl);
    
    typedef Parameters<real> Par_t;
    typedef EvolveSurface<real, DT, the_device> Evolve_t;
    typedef Evolve_t::Vec_t Vec_t;

    // Setting the parameters
    Par_t sim_par;
    sim_par.n_surfs = 2;   
    sim_par.ts = 1;    
    sim_par.time_horizon = 2;
    sim_par.scheme = Explicit;
    sim_par.bg_flow_param = 0.1;
    sim_par.rep_maxit = 20;
    sim_par.save_data = true;    
    sim_par.save_stride = 1;
    sim_par.save_file_name = "TaylorVortex.out";
    COUT(sim_par);
    
    //Cleaning the slate
    remove(sim_par.save_file_name.c_str());


    //Initial vesicle positions 
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
    DataIO myIO;
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0, 0, x0.getSubLength());
    
    //Making Centers And Populating The Prototype
    int nVec = sim_par.n_surfs;
    real* cntrs_host =  new real[nVec * DIM];
    for(int ii=0; ii<nVec; ++ii)
    {
        cntrs_host[DIM*ii    ] = 0;
        cntrs_host[DIM*ii + 1] = 3*ii;
        cntrs_host[DIM*ii + 2] = 0;
    }
    Array<real, DT, the_device> cntrs(DIM * sim_par.n_surfs);
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    Populate(x0, cntrs);

    //Reading Operators From File
    bool readFromFile = true;
    Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    TaylorVortex<Vec_t> vInf(sim_par.bg_flow_param);
   
    //Finally, Evolve surface
    PROFILESTART();    
    Evolve_t Es(sim_par, Mats, x0, &vInf, &StokesAlltoAll);
    Es.Evolve();    
    PROFILEEND("",0);
    PROFILEREPORT(SortFlop);
}



