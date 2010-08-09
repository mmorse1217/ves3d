#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
#include "EvolveSurface.h"

//CPU
extern const Device<CPU> dev0(0);

//GPU
extern const Device<GPU> dev1(0);
extern const Device<GPU> dev2(1);

typedef float real;
typedef VesInteraction<real> Interaction_t;
typedef RepartitionGateway<real> Repart_t;

template<enum DeviceType DT, const Device<DT> &dev>
void setupES(Parameters<real> &sim_par, 
    Interaction_t& interaction, Repart_t &repart)
{
    typedef Scalars<real, DT, dev> Sca_t;
    typedef Vectors<real, DT, dev> Vec_t;
    typedef Surface<Sca_t,Vec_t> Sur_t;
    
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //IO
    DataIO myIO;

    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0, 0, x0.getSubLength());
    
    
    //Making Centers And Populating The Prototype
    int nVec = sim_par.n_surfs;
    real* cntrs_host =  new real[nVec * DIM];
    for(int ii=0; ii<nVec; ++ii)
    {
        cntrs_host[DIM*ii    ] = 5 * omp_get_thread_num();
        cntrs_host[DIM*ii + 1] = 5 * ii;
        cntrs_host[DIM*ii + 2] = 0;
    }
    
    Array<real, DT, dev> cntrs(DIM * sim_par.n_surfs);
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    Populate(x0, cntrs);
    delete[] cntrs_host;

    //Reading Operators From File
    bool readFromFile = true;
    OperatorsMats<Sca_t> Mats(readFromFile, sim_par);
    
    //Making The Surface, And Time Stepper
    Sur_t S(x0, Mats);
    Monitor<Sur_t> M(sim_par);
    EvolveSurface<Sur_t, Interaction_t, Monitor<Sur_t>, Repart_t> Es(Mats, sim_par, M, repart);
    Es(S, interaction);
}
    
int main(int argc, char **argv)
{
    typedef Parameters<real> Par_t;
    
    int nThreads = 3;
    Interaction_t interaction(&StokesAlltoAll, nThreads);
    Repart_t repart(NULL, nThreads);
    
    // Making multiple threads
#pragma omp parallel num_threads(nThreads)
    {
        // Setting the parameters
        Par_t sim_par;
        sim_par.n_surfs = 128; 
        sim_par.ts = .1;    
        sim_par.time_horizon = 3;
        sim_par.scheme = Explicit;
        sim_par.bg_flow_param = 0.1;
        sim_par.rep_maxit = 20;
        sim_par.save_data = false;
        sim_par.save_stride = 1;
        sim_par.save_file_name = "";
        //COUT(sim_par);
                
        DataIO myIO;
      
        switch ( omp_get_thread_num() )
        {
            case 0:
                setupES<CPU, dev0>(sim_par, interaction, repart);
                break;

            case 1:
                setupES<GPU, dev1>(sim_par, interaction, repart);
                break;

            case 2:
                setupES<GPU, dev2>(sim_par, interaction, repart);
                break;
        }
    }
}

