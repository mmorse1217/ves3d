#include "EvolveSurface.h"
#include "Device.h"
#include "Scalars.h"
#include "Vectors.h"
#include "SurfaceParams.h"
#include "OperatorsMats.h"
#include "Surface.h"

extern const Device<CPU> the_cpu_dev(0);

typedef double real;

int main(int argc, char **argv)
{
    typedef containers::Scalars<real, CPU, the_cpu_dev> Sc;
    typedef containers::Vectors<real, CPU, the_cpu_dev> Vc;
    typedef Surface<Sc,Vc> Sur;

    int p(12), nSur(1);
    //int fLen(2*p*(p+1));
    int dLen(6*p*(p+1));
    
    SurfaceParams<Sc::value_type> par; 
    par.p_ = p;
    par.n_surfs_ = 1;
    par.kappa_ = 1e-2;
    par.filter_freq_ = 2*p/3;
    par.rep_ts_ = 1e-1;
    par.rep_max_vel_ = 1e-1;
    par.rep_iter_max_ = 100;
    par.rep_up_freq_ = 2*p;
    par.rep_filter_freq_ = p/3;
    
    //IO
    DataIO<Sc::value_type, CPU> myIO(the_cpu_dev);
    
    // initializing vesicle positions from text file
    Vc x0(nSur, p);
                       
    char fname[400];
    sprintf(fname,"%s/precomputed/dumbbell_cart12_single.txt",
        getenv("VES3D_DIR"));
    myIO.ReadData(fname,x0.size(),x0.begin());
    
    Sur S(x0);  

    real dt = .01;
    real T = 2 * dt;    
    EvolveSurface<Sur> Es;
    Es(S, T, dt);

    myIO.WriteData("EvolveSurf.txt", S.getPosition().size(), S.getPosition().begin());
}
