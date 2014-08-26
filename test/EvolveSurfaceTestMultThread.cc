#include "EvolveSurface.h"

extern const Device<CPU> device(0);

typedef float real;
typedef EvolveSurface<real, CPU, device> Evolve_t;
typedef typename Evolve_t::Vec_t Vec_t;
typedef typename Evolve_t::Interaction_t Interaction_t;
typedef typename Evolve_t::Repartition_t Repartition_t;

#ifndef Doxygen_skip
void MyRepart(size_t nv, size_t stride, const real* x,
    const real* tension, size_t* nvr,
    real** xr, real** tensionr, void*)
{
    int expRate(2);
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

template<enum DeviceType DT, const Device<DT> &DEVICE>
void EvolveSurfaceTestMT(Parameters<real> &sim_par, int nsurfs,
    Interaction_t* interaction, Repartition_t* repartition)
{
    //Initial vesicle positions
    Vec_t x0(nsurfs, sim_par.sh_order);

    //reading the prototype form file
    DataIO myIO;
    char fname[300];
    string prec = (typeid(real) == typeid(float)) ? "float" : "double";
    sprintf(fname,"precomputed/dumbbell_%u_%s.txt",sim_par.sh_order,prec.c_str());
    myIO.ReadData(fname, x0, 0, x0.getSubLength());

    //Making Centers And Populating The Prototype there is no
    assert(interaction == NULL);
    Array<real, DT, DEVICE> cntrs(DIM * nsurfs);
    {
        real* cntrs_host =  new real[cntrs.size()];
        for(int ii=0; ii<cntrs.size(); ++ii)
            cntrs_host[ii] = 0;

        cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
            cntrs.size() * sizeof(real), MemcpyHostToDevice);
        delete[] cntrs_host;
    }
    Populate(x0, cntrs);

    //Reading Operators From File
    bool readFromFile = true;
    typename Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

    //Finally, Evolve surface
    Evolve_t Es(sim_par, Mats, x0, &vInf, NULL, interaction, repartition);

    QC ( Es.Evolve() );

    //checking that everything is still identical
    assert(interaction == NULL);
    cntrs.resize(DIM * Es.S_->getNumberOfSurfaces());
    x0.replicate(Es.S_->getPosition());
    axpy(0, x0, Es.S_->getPosition(), x0);
    {
        real* cntrs_host =  new real[cntrs.size()];
        for(int ii=0; ii<cntrs.size(); ++ii)
            cntrs_host[ii] = 0;

        cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
            cntrs.size() * sizeof(real), MemcpyHostToDevice);
        delete[] cntrs_host;
    }
    Populate(x0, cntrs);
    axpy(-1,x0, Es.S_->getPosition(), x0);

#pragma omp barrier
#pragma omp critical
    COUT("\n\n  Thread "<<omp_get_thread_num()<<" error : "<<MaxAbs(x0)<<endl);

#pragma omp barrier
    return;
}
#endif Doxygen_skip

int main(int argc, char **argv)
{
    COUT("\n\n =========================================\n"
        <<"  EvolveSurface with repartitioning test: "
        <<"\n ========================================="<<endl);

    int nThreads = 3;
    Repartition_t repartition(&MyRepart);
    // Making multiple threads
#pragma omp parallel num_threads(nThreads)
    {
        // Setting the parameters
        Parameters<real> sim_par;

        sim_par.n_surfs = 0;
        sim_par.sh_order = 12;
        sim_par.filter_freq = 8;
        sim_par.rep_up_freq = 24;
        sim_par.rep_filter_freq = 4;

        sim_par.ts = 1;
        sim_par.time_horizon = 3;
        sim_par.scheme = BlockImplicit; //Explicit
        sim_par.singular_stokes = Direct;
        sim_par.bg_flow_param = 0;
        sim_par.rep_maxit = 20;

        int nsurfs(10*(omp_get_thread_num()+1));
        EvolveSurfaceTestMT<CPU, device>(sim_par, nsurfs, NULL, &repartition);
    }
}
