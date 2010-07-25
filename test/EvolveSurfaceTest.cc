#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
#include "OperatorsMats.h"
#include "EvolveSurface.h"

extern const Device<CPU> the_cpu_device(0);
extern const Device<GPU> the_gpu_device(0);

typedef double real;
typedef double fmm_value_type;

#ifndef Doxygen_skip

template<enum DeviceType DT, const Device<DT> &DEVICE>
void EvolveSurfaceTest(Parameters<real> &sim_par)
{
    typedef Scalars<real, CPU, DEVICE> Sca_t;
    typedef Vectors<real, CPU, DEVICE> Vec_t;

    typedef Surface<Sca_t,Vec_t> Sur_t;
    typedef VesInteraction<fmm_value_type> Interaction_t;
    
    //IO
    DataIO myIO;

    //Initializing vesicle positions from text file
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
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
    Array<real, DT, DEVICE> cntrs(DIM * sim_par.n_surfs);
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    Populate(x0, cntrs);

    // The Interaction Class
    Interaction_t Interaction(&StokesAlltoAll);

    //Reading Operators From File
    bool readFromFile = true;
    OperatorsMats<Sca_t> Mats(readFromFile, sim_par);

    //Making The Surface, And Time Stepper
    Sur_t S(x0, Mats);
    Monitor<Sur_t> M(sim_par);
    RepartitionGateway<real> repart(NULL);
    EvolveSurface<Sur_t, Interaction_t> Es(Mats, sim_par, M, repart);
   
    Es(S, Interaction);
}
#endif //Doxygen_skip

int main(int argc, char **argv)
{
    COUT("\n\n ========================\n  EvolveSurface test: "
        <<"\n ========================"<<endl);

    typedef Parameters<real> Par_t;
    // Setting the parameters
    Par_t sim_par;
    sim_par.n_surfs = 0;   
    sim_par.ts = 1;    
    sim_par.time_horizon = 5;
    sim_par.scheme = Explicit;
    sim_par.bg_flow_param = 0.1;
    sim_par.rep_maxit = 20;
    sim_par.save_data = true;    
    sim_par.save_stride = 1;
    sim_par.save_file_name = "EvolveSurf.out";
    COUT(sim_par);
    
    //Cleaning the slate
    remove(sim_par.save_file_name.c_str());

    COUT("\n ------------ \n  CPU device: \n ------------"<<endl);
    EvolveSurfaceTest<CPU,the_cpu_device>(sim_par);
    PROFILEREPORT(SortTime);    

#ifdef GPU_ACTIVE
    PROFILECLEAR();
    COUT("\n ------------ \n  GPU device: \n ------------"<<endl);
    EvolveSurfaceTest<GPU, the_gpu_device>(sim_par);
    PROFILEREPORT(SortTime);

#endif //GPU_ACTIVE


}


