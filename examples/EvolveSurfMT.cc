#include "EvolveSurface.h"

extern const Device<CPU> dev0(0);
extern const Device<GPU> dev1(0);

typedef float real;
typedef VesInteraction<real> Interaction_t;
typedef Repartition<real> Repart_t;

void MyRepart(size_t nv, size_t stride, const real* x, 
    const real* tension, size_t* nvr, 
    real** xr, real** tensionr, void*)
{
    int expRate(1);
    *nvr = expRate*nv;
    *xr = new real[*nvr * stride * DIM];
    *tensionr = new real[*nvr * stride];
        
    int ll(nv * stride);
    for(int ii=0; ii<expRate; ++ii)
    {
        memcpy(*xr       + ii * ll * DIM, x       , DIM * ll * sizeof(real));
        memcpy(*tensionr + ii * ll      , tension ,       ll * sizeof(real));
    }
}

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
    char fname[300];
    string prec = (typeid(real) == typeid(float)) ? "float" : "double"; 
    sprintf(fname,"precomputed/dumbbell_%u_%s.txt",sim_par.sh_order,prec.c_str());
    myIO.ReadData(fname, x0, 0, x0.getSubLength());
    
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
    

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

    //Finally, Evolve surface
    EvolveSurface<real, DT, dev, Interaction_t, Repart_t> Es(sim_par, Mats, 
        x0, &vInf, NULL, &interaction, &repart);
    
    QC ( Es.Evolve() );
}
    
int main(int argc, char **argv)
{
    typedef Parameters<real> Par_t;
    
    int nThreads = 2;
    Interaction_t interaction(&StokesAlltoAll);
    Repart_t repart(&MyRepart);
    
    // Making multiple threads
#pragma omp parallel num_threads(nThreads)
    {
        // Setting the parameters
        Par_t sim_par;
        sim_par.n_surfs = 10; 
        sim_par.ts = .1;    
        sim_par.time_horizon = .3;
        sim_par.scheme = Explicit;
        sim_par.bg_flow_param = 0.1;
        sim_par.rep_maxit = 20;
        sim_par.save_data = false;
        sim_par.save_stride = 1;
        sim_par.save_file_name = "";
        sim_par.singular_stokes = Direct;
                
        DataIO myIO;
      
        switch ( omp_get_thread_num() )
        {
            case 0:
                setupES<CPU, dev0>(sim_par, interaction, repart);
                break;

            case 1:
                setupES<GPU, dev1>(sim_par, interaction, repart);
                break;
        }
    }
}

