#include "EvolveSurface.h"
#include "CPUKernels.h"

extern const Device<CPU> the_cpu_device(0);
extern const Device<GPU> the_gpu_device(0);

typedef float real;

template<enum DeviceType DT, const Device<DT> &DEVICE>
void EvolveSurfaceTest(Parameters<real> &sim_par)
{
    typedef EvolveSurface<real, DT, DEVICE> Evolve_t;
    typedef typename Evolve_t::Vec_t Vec_t;

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
    Array<real, DT, DEVICE> cntrs(DIM * sim_par.n_surfs);
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    Populate(x0, cntrs);

    //Reading Operators From File
    bool readFromFile = true;
    typename Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

    //Finally, Evolve surface
    Evolve_t Es(sim_par, Mats, x0, &vInf);//&StokesAlltoAll);
    Es.Evolve();
}

int main(int argc, char **argv)
{
    COUT("\n\n ========================\n  EvolveSurface test: "
        <<"\n ========================"<<endl);

    typedef Parameters<real> Par_t;
    // Setting the parameters
    Par_t sim_par;
    
    sim_par.sh_order = 24;
    sim_par.filter_freq = 16;
    sim_par.rep_up_freq = 48;
    sim_par.rep_filter_freq = 8;
 
    sim_par.n_surfs = 1024;   
    sim_par.ts = 1;    
    sim_par.time_horizon = 1;
    sim_par.scheme = Explicit;
    sim_par.bg_flow_param = 0.1;
    sim_par.rep_maxit = 20;
    sim_par.save_data = false;    
    sim_par.save_stride = 1;
    sim_par.save_file_name = "EvolveSurf.out";
    COUT(sim_par);
    
    //Cleaning the slate
    remove(sim_par.save_file_name.c_str());
    
//     PROFILESTART();
//     COUT("\n ------------ \n  CPU device: \n ------------"<<endl);
//     EvolveSurfaceTest<CPU,the_cpu_device>(sim_par);
//     PROFILEEND("",0);
//     PROFILEREPORT(SortFlopRate);    

#ifdef GPU_ACTIVE
    PROFILECLEAR();
    COUT("\n ------------ \n  GPU device: \n ------------"<<endl);
    EvolveSurfaceTest<GPU, the_gpu_device>(sim_par);
    PROFILEREPORT(SortFlopRate);

#endif //GPU_ACTIVE
}


