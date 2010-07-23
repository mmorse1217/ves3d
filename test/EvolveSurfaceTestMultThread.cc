#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
#include "EvolveSurface.h"

extern const Device<CPU> dev1(0);
extern const Device<CPU> dev2(0);

typedef double real;
typedef VesInteraction<real> Interaction_t;
typedef RepartitionGateway<real> Repart_t;

#ifndef Doxygen_skip

void MyRepart(size_t nv, size_t stride, const double* x, 
    const double* tension, size_t* nvr, 
    double** xr, double** tensionr, void*)
{
    *nvr = nv;
    *xr = new double[nv * stride * DIM];
    *tensionr = new double[nv * stride];
    
    memcpy(*xr      ,       x, DIM * nv * stride * sizeof(double));
    memcpy(*tensionr, tension,       nv * stride * sizeof(double));
}

template<enum DeviceType DT, const Device<DT> &DEVICE>
void EvolveSurfaceTest(Parameters<real> &sim_par, 
    Interaction_t& Interaction, Repart_t &repart)
{
    typedef Scalars<real, CPU, DEVICE> Sca_t;
    typedef Vectors<real, CPU, DEVICE> Vec_t;

    typedef Surface<Sca_t,Vec_t> Sur_t;
    
    //IO
    DataIO myIO;

    //Initializing vesicle positions from text file
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",
        x0, 0, DIM*x0.getStride());
    
    //Making Centers And Populating The Prototype
    int nVec = sim_par.n_surfs;
    real* cntrs_host =  new real[nVec * DIM];
    for(int ii=0; ii<nVec; ++ii)
    {
        cntrs_host[DIM*ii    ] = 3 * omp_get_thread_num();
        cntrs_host[DIM*ii + 1] = 3*ii;
        cntrs_host[DIM*ii + 2] = 0;
    }
    Array<real, DT, DEVICE> cntrs(DIM * sim_par.n_surfs);
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    Populate(x0, cntrs);

    //Reading Operators From File
    bool readFromFile = true;
    OperatorsMats<Sca_t> Mats(readFromFile, sim_par);

    //Making The Surface, And Time Stepper
    Sur_t S(x0, Mats);
    Monitor<Sur_t> M(sim_par);
    EvolveSurface<Sur_t, Interaction_t, Monitor<Sur_t>, Repart_t> Es(Mats, sim_par, M, repart);
    Es(S, Interaction);
    
}
#endif //Doxygen_skip

int main(int argc, char **argv)
{
    typedef Parameters<real> Par_t;

    int nThreads = 2;
    Interaction_t interaction(&StokesAlltoAll, nThreads);
    Repart_t repart(&MyRepart, nThreads);
    
    // Making multiple threads
#pragma omp parallel num_threads(nThreads)
    {
        // Setting the parameters
        Par_t sim_par;
        sim_par.ts = 1;    
        sim_par.time_horizon = 3;
        sim_par.scheme = Explicit;
        sim_par.bg_flow_param = 0.1;
        sim_par.rep_maxit = 20;
        sim_par.save_data = true;    
        sim_par.save_stride = 1;
        COUT(sim_par);
        remove(sim_par.save_file_name.c_str());
        
        if(omp_get_thread_num() == 0)
        {   
            sim_par.save_file_name = "EvolveSurfMT_0.out";
            sim_par.n_surfs = 1; 
            EvolveSurfaceTest<CPU,dev1>(sim_par, interaction, repart);
        }
        else
        {
            sim_par.save_file_name = "EvolveSurfMT_1.out";
            sim_par.n_surfs = 3;
            EvolveSurfaceTest<CPU,dev2>(sim_par, interaction, repart);
        }        
    }
    
    PROFILEREPORT(SortTime);
}
