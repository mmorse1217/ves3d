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

template<typename Sca, typename Vec, enum DeviceType DT>
void EvolveSurfTest(const Device<DT> &dev, Parameters<real> &sim_par)
{
    typedef Surface<Sca,Vec> Sur_t;
    typedef VesInteraction<fmm_value_type> Interaction_t;
    
    //IO
    DataIO myIO;

    //Initializing vesicle positions from text file
    Vec x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0, 0, DIM*x0.getStride());
    
    //Making Centers And Populating The Prototype
    if ( sim_par.n_surfs > 1 )
    {
        Vec cntrs(x0.getNumSubs(), 0, make_pair(1,1));
        real cntrs_host[] = {-5, 0,  1,
                             5, 0, -1};
        
        cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
            cntrs.size() * sizeof(real), MemcpyHostToDevice);
        
        Populate(x0, cntrs);
    }

    // The Interaction Class
    Interaction_t Interaction(NULL);//&StokesAlltoAll);

    //Reading Operators From File
    bool readFromFile = true;
    OperatorsMats<Sca> Mats(readFromFile, sim_par);

    //Making The Surface, And Time Stepper
    Sur_t S(x0, Mats);
    Monitor<Sur_t> M(sim_par);
    EvolveSurface<Sur_t, Interaction_t> Es(Mats, sim_par, M);
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
    sim_par.n_surfs = 1;   
    sim_par.ts = 1;    
    sim_par.time_horizon = 5;
    sim_par.scheme = Explicit;//SemiImplicit;
    sim_par.bg_flow_param = 0.1;
    sim_par.rep_maxit = 20;
    sim_par.save_data = true;    
    sim_par.save_stride = 1;
    sim_par.save_file_name = "EvolveSurf.out";
    COUT(sim_par);
    
    //Cleaning the slate
    remove(sim_par.save_file_name.c_str());

    COUT("\n ------------ \n  CPU device: \n ------------"<<endl);

    typedef Scalars<real, CPU, the_cpu_device> ScaCPU_t;
    typedef Vectors<real, CPU, the_cpu_device> VecCPU_t;
    
    EvolveSurfTest<ScaCPU_t, VecCPU_t, CPU>(the_cpu_device, sim_par);
    PROFILEREPORT(SortTime);    

#ifdef GPU_ACTIVE
    PROFILECLEAR();
    COUT("\n ------------ \n  GPU device: \n ------------"<<endl);
    typedef Scalars<real, GPU, the_gpu_device> ScaGPU_t;
    typedef Vectors<real, GPU, the_gpu_device> VecGPU_t;
     
    EvolveSurfTest<ScaGPU_t, VecGPU_t, GPU>(the_gpu_device, sim_par);

    PROFILEREPORT(SortTime);
#endif //GPU_ACTIVE


}


