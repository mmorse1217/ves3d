#include "EvolveSurface.h"
#include "Error.h"
#include "HelperFuns.h"
#include "CPUKernels.h"

typedef double real;
typedef Device<CPU> DevCPU;
typedef Device<GPU> DevGPU;

extern const DevCPU the_cpu_device(0);
extern const DevGPU the_gpu_device(0);

template<typename DT, const DT &DEVICE>
void EvolveSurfaceTest(Parameters<real> &sim_par)
{
    typedef EvolveSurface<real, DT, DEVICE> Evolve_t;
    typedef typename Evolve_t::Vec_t Vec_t;

    //Initial vesicle positions
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);

    //reading the prototype form file
    DataIO myIO(sim_par.checkpoint_file_name,DataIO::ASCII);
    char fname[300];
    std::string prec = (typeid(real) == typeid(float)) ? "float" : "double";
    sprintf(fname,"precomputed/dumbbell_%u_%s.txt",sim_par.sh_order,prec.c_str());
    myIO.ReadData(FullPath(fname), x0, DataIO::ASCII, 0, x0.getSubLength());

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
        cntrs.size() * sizeof(real), DT::MemcpyHostToDevice);
    Populate(x0, cntrs);

    //Reading Operators From File
    bool readFromFile = true;
    typename Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    BgFlowBase<Vec_t> *vInf(NULL);
    CHK(BgFlowFactory(sim_par, &vInf));

    typename Evolve_t::Interaction_t interaction(&StokesAlltoAll);

    //Finally, Evolve surface
    Evolve_t Es(&sim_par, Mats, vInf, NULL, &interaction, NULL, NULL, &x0);

    CHK ( Es.Evolve() );

    //testing streaming
    std::stringstream s1,s2;
    Es.pack(s1, Streamable::ASCII);
    Evolve_t E2(&sim_par, Mats, vInf, NULL, &interaction);
    E2.unpack(s1, Streamable::ASCII);
    E2.pack(s2, Streamable::ASCII);
    ASSERT(s1.str()==s2.str(), "bad streaming");
}

int main(int argc, char **argv)
{
    COUT("========================\n  EvolveSurface test: "
        <<"\n========================");

    typedef Parameters<real> Par_t;
    // Setting the parameters
    Par_t sim_par;

    sim_par.sh_order        = 6;
    sim_par.filter_freq     = 4;
    sim_par.rep_up_freq     = 12;
    sim_par.rep_filter_freq = 2;

    sim_par.n_surfs              = 2;
    sim_par.ts                   = 1;
    sim_par.time_horizon         = 2;
    sim_par.scheme               = JacobiBlockExplicit;
    sim_par.singular_stokes      = Direct;
    sim_par.bg_flow_param        = 0;
    sim_par.upsample_interaction = true;
    sim_par.rep_maxit            = 20;
    sim_par.checkpoint           = true;
    sim_par.checkpoint_stride    = 1;
    sim_par.checkpoint_file_name = "EvolveSurf.out";
    COUT(sim_par);

    //Cleaning the slate
    remove(sim_par.checkpoint_file_name.c_str());

    CLEARERRORHIST();
    PROFILESTART();
    COUTDEBUG("Testing CPU device");
    EvolveSurfaceTest<DevCPU,the_cpu_device>(sim_par);
    PROFILEEND("",0);
    PRINTERRORLOG();
    PROFILEREPORT(SortFlopRate);

#ifdef GPU_ACTIVE
    CLEARERRORHIST();
    PROFILECLEAR();
    COUTDEBUG("Testing GPU device");
    EvolveSurfaceTest<DevGPU, the_gpu_device>(sim_par);
    PRINTERRORLOG();
    PROFILEREPORT(SortFlopRate);

#endif //GPU_ACTIVE
}
