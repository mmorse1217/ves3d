#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
#include "OperatorsMats.h"
#include "EvolveSurface.h"

typedef float real;
#define DT GPU
extern const Device<DT> the_device(0);

template<const Device<DT> &DEVICE>
void EvolveSurfaceTest(Parameters<real> &sim_par)
{
    typedef Scalars<real, DT, DEVICE> Sca_t;
    typedef Vectors<real, DT, DEVICE> Vec_t;

    typedef Surface<Sca_t,Vec_t> Sur_t;
    typedef VesInteraction<real> Interaction_t;
    
    //IO
    PROFILESTART();
    DataIO myIO;

    //Initializing vesicle positions from text file
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
    char fname[300];
    sprintf(fname, "precomputed/biconcave_ra85_%u",
        sim_par.sh_order);
    myIO.ReadData(fname, x0, 0, x0.getSubLength());
    
    //Making Centers And Populating The Prototype
    int nVec = sim_par.n_surfs;
    real* cntrs_host =  new real[nVec * DIM];
    for(int ii=0; ii<nVec; ++ii)
    {
        cntrs_host[DIM*ii    ] = 0;
        cntrs_host[DIM*ii + 1] = 5*ii;
        cntrs_host[DIM*ii + 2] = 0;
    }
    Array<real, DT, DEVICE> cntrs(DIM * sim_par.n_surfs);
    cntrs.getDevice().Memcpy(cntrs.begin(), cntrs_host,
        cntrs.size() * sizeof(real), MemcpyHostToDevice);
    Populate(x0, cntrs);

    // The Interaction Class
    Interaction_t Interaction(NULL);

    //Reading Operators From File
    bool readFromFile = true;
    OperatorsMats<Sca_t> Mats(readFromFile, sim_par);

    //Making The Surface, And Time Stepper
    Sur_t S(x0, Mats);
    Monitor<Sur_t> M(sim_par);
    RepartitionGateway<real> repart(NULL);
    EvolveSurface<Sur_t, Interaction_t> Es(Mats, sim_par, M, repart);
 
    PROFILEEND("setup_",0);

    Es(S, Interaction);
}

int main(int argc, char **argv)
{
    typedef Parameters<real> Par_t;
 
    // Setting the parameters
    Par_t sim_par;
    sim_par.n_surfs = 0;   
    sim_par.ts = .1;    
    sim_par.time_horizon = .3;
    sim_par.bg_flow_param = 0.1;
    sim_par.rep_maxit = 20;
    sim_par.save_data = false;    
    
    sim_par.scheme = Explicit;
    int maxexp = 11; // nmax = 4096
    int p[] = {6, 12, 16, 24};
    int plength = 5;

    for(int ii=0;ii<plength; ++ii)
    {
        int n0 = 1;
        sim_par.sh_order = p[ii];
        sim_par.filter_freq = 2*p[ii]/3;
        sim_par.rep_up_freq = 2*p[ii];
        sim_par.rep_filter_freq = p[ii]/3;
        //rep_ts(1),
        for(int jj = 0; jj<maxexp; ++jj)
        {
            PROFILECLEAR();
            PROFILESTART();
            n0 *= 2;
            sim_par.n_surfs = n0; 
            
            COUT("\n --- n = " << n0 <<" p = "
                <<p[ii]<<"------------------------------------------------------------"<<endl);
            COUT(sim_par);
            
            EvolveSurfaceTest<the_device>(sim_par);
            PROFILEEND("",0);
            PROFILEREPORT(SortTime);
        }
    }
}


