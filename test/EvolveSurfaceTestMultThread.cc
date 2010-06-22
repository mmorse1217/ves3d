#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
//#include "EvolveSurface.h"

extern const Device<CPU> the_cpu_dev(0);
extern const Device<GPU> the_gpu_dev(0);

typedef float real;
typedef double fmm_value_type;

void FMM(const fmm_value_type*src, const fmm_value_type*den, size_t np, fmm_value_type*pot)
{
    cout<<"FMM with "<<np<<endl;
    sleep(2);
}

int main(int argc, char **argv)
{
    typedef containers::Scalars<real, CPU, the_cpu_dev> ScaCPU;
    typedef containers::Vectors<real, CPU, the_cpu_dev> VecCPU;
    typedef Surface<ScaCPU,VecCPU> SurCPU;

    typedef containers::Scalars<real, GPU, the_gpu_dev> ScaGPU;
    typedef containers::Vectors<real, GPU, the_gpu_dev> VecGPU;
    typedef Surface<ScaCPU,VecCPU> SurGPU;

    typedef Parameters<real> Par;

    // Setting the parameters
    Par::getInstanceModifiable().n_surfs = 2;   
    //Par::getInstanceModifiable().n_steps = 20;
    Par::getInstanceModifiable().ts = .5;    
    Par::getInstanceModifiable().time_horizon = 30;
    Par::getInstanceModifiable().inner_solver_maxit = 15;    
    //Par::getInstanceModifiable().bg_flow_param = 0;    
    cout<<Par::getInstance()<<endl;

    //IO
    DataIO<real, CPU> myIO(the_cpu_dev);
    
    //Initializing vesicle positions from text file
    Vec x0(Par::getInstance().n_surfs,
        Par::getInstance().sh_order);
    
    //reading the prototype form file
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0.getSubN(1)-x0.begin(), x0.begin());
    
    //Making centers and populating the prototype
    Vec cntrs(x0.getNumSubs(), 0, make_pair(1,1));
    real cntrs_host[] = {-5, 0,  1,
                         5, 0, -1};
    
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    
    Populate(x0, cntrs);

    // The interaction class
    int nThreads(2);
    VesInteraction<fmm_value_type> interaction(nThreads, &FMM);

// Making multiple threads
#pragma omp parallel num_threads(nThreads)
    {
        if(omp_get_thread_num() == 0)
        {   
            VecCPU thVec(omp_get_thread_num() + 1, 4);
            interaction(thVec, thVec, thVec);
            cout<<"Thread "<<omp_get_thread_num()<<endl;
        }
        else
        {
            VecGPU thVec(omp_get_thread_num() + 2, 4);
            interaction(thVec, thVec, thVec);
            cout<<"Thread "<<omp_get_thread_num()<<endl;
        }        
    }



//     //Making the surface, and time stepper
//     Sur S(x0);
//     EvolveSurface<Sur> Es;
//     Es(S);

//     myIO.WriteData("EvolveSurf.txt", S.getPosition().size(), 
//         S.getPosition().begin());

//     PROFILEREPORT(SortTime);
}
