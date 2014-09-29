#undef PROFILING // profiling isn't for multithreaded (we're imitating mpi here)

#include "EvolveSurface.h"
#include "Error.h"

typedef double real;
typedef Device<CPU> DevCPU;
typedef Device<GPU> DevGPU;

extern const DevCPU the_cpu(0);
extern const DevGPU the_gpu(0);

#ifndef Doxygen_skip
void MyRepart(size_t nv, size_t stride, const real* x,
    const real* tension, size_t* nvr,
    real** xr, real** tensionr, void*)
{
    int expRate(2);
    *nvr = expRate*nv;
    INFO("Increasing number of surfaces from "<<nv<<" to "<<*nvr);
    *xr = new real[*nvr * stride * DIM];
    *tensionr = new real[*nvr * stride];

    int ll(nv * stride);
    for(int ii=0; ii<expRate; ++ii)
    {
        memcpy(*xr       + ii * ll * DIM, x       , DIM * ll * sizeof(real));
        memcpy(*tensionr + ii * ll      , tension ,       ll * sizeof(real));
    }
}

template<typename T,typename Evolve_t>
void EvolveEachObj(Parameters<T> &sim_par,
    int nsurfs,
    typename Evolve_t::Interaction_t* interaction,
    typename Evolve_t::Repartition_t* repartition)
{
    ASSERT(interaction == NULL, "Doesn't support interaction yet");

    typedef typename Evolve_t::Arr_t Arr_t;
    //    typedef typename Evolve_t::Sca_t Sca_t;
    typedef typename Evolve_t::Vec_t Vec_t;
    typedef typename Evolve_t::device_type DT;
    typedef typename Evolve_t::Mats_t Mats_t;
    //    typedef typename Evolve_t::Interaction_t Interaction_t;
    //typedef typename Evolve_t::Repartition_t Repartition_t;

    //Initial vesicle positions
    Vec_t x0(nsurfs, sim_par.sh_order);
    DataIO myIO;
    char fname[300];
    std::string prec = (typeid(T) == typeid(float)) ? "float" : "double";
    sprintf(fname,"%s/precomputed/dumbbell_%u_%s.txt",VES3D_PATH,sim_par.sh_order,prec.c_str());

    //reading the prototype form file (critical only for the log file)
#pragma omp critical
    myIO.ReadData(fname, x0, DataIO::ASCII, 0, x0.getSubLength());

    //Making Centers And Populating The Prototype
    int nVec = nsurfs;
    Arr_t cntrs(DIM * nsurfs);
    T* cntrs_host =  new T[nVec * DIM];

    for(int ii=0; ii<cntrs.size(); ++ii)
        cntrs_host[ii] = 0;

    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), DT::MemcpyHostToDevice);
    Populate(x0, cntrs);

    //Reading Operators From File
    bool readFromFile = true;
    std::auto_ptr<Mats_t> Mats;

#pragma omp critical  //just to get a nicer log
    Mats.reset(new Mats_t(readFromFile, sim_par));

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

    // Evolver
    Evolve_t Es(sim_par, *Mats, x0, &vInf, NULL, interaction, repartition);

    CHK ( Es.Evolve() );

#pragma omp barrier
#pragma omp critical //just to get a nicer log
    {
        INFO("Checking thread "<<omp_get_thread_num());
        //checking that everything is  identical
        ASSERT(Es.S_->getNumberOfSurfaces()==80,
            "Uniform distribution : " <<Es.S_->getNumberOfSurfaces());

        cntrs.resize(DIM * Es.S_->getNumberOfSurfaces());
        x0.replicate(Es.S_->getPosition());
        axpy(0, x0, Es.S_->getPosition(), x0);
        {
            real* cntrs_host =  new real[cntrs.size()];
            for(int ii=0; ii<cntrs.size(); ++ii)
                cntrs_host[ii] = 0;

            cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
                cntrs.size() * sizeof(real), DT::MemcpyHostToDevice);
            delete[] cntrs_host;
        }
        Populate(x0, cntrs);
        axpy(-1,x0, Es.S_->getPosition(), x0);

        ASSERT(MaxAbs(x0)<1e-14,"All surfaces are the same");

        INFO(emph<<"*** Thread "<<omp_get_thread_num()
            <<" passed, error : "<<MaxAbs(x0)<<" ***"<<emph);
    }
}

template<typename T, typename DT, const DT &device>
void TestEvolveSurfMT(){

    INFO("Testing "<<DT::type()<<" device");

    typedef EvolveSurface<T, DT, device> Evolve_t;
    typedef typename Evolve_t::Interaction_t Interaction_t;
    typedef typename Evolve_t::Repartition_t Repartition_t;

    Repartition_t repartition(&MyRepart);

    // Setting the parameters
    int nThreads = 3;
    Parameters<T> sim_par;

    sim_par.sh_order        = 6;
    sim_par.filter_freq     = 4;
    sim_par.rep_up_freq     = 12;
    sim_par.rep_filter_freq = 2;

    sim_par.ts                   = 1;
    sim_par.time_horizon         = 2;
    sim_par.scheme               = BlockImplicit; //Explicit
    sim_par.singular_stokes      = Direct;
    sim_par.bg_flow_param        = 0;
    sim_par.upsample_interaction = true;
    sim_par.rep_maxit            = 20;
    sim_par.save_data            = false;
    COUT(sim_par);

    // Making multiple threads
#pragma omp parallel num_threads(nThreads)
    {
        int nsurfs = 10*(omp_get_thread_num()+1);
        EvolveEachObj<T,Evolve_t>(sim_par, nsurfs, NULL, &repartition);
    }
}
#endif Doxygen_skip

int main(int argc, char **argv)
{
    COUT("=========================================\n"
        <<" EvolveSurface with repartitioning test:\n"
        <<"=========================================");

    TestEvolveSurfMT<real,DevCPU,the_cpu>();
#ifdef GPU_ACTIVE
    TestEvolveSurfMT<real,DevGPU,the_gpu>();
#endif //GPU_ACTIVE
}
