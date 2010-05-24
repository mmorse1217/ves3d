#include "Device.h"
#include "Scalars.h"
#include "Vectors.h"
#include "SurfaceParams.h"
#include "OperatorsMats.h"
#include "Surface.h"
#include "CurvatureFlow.h"

extern const Device<CPU> the_cpu_dev(0);

int main(int argc, char **argv)
{
    typedef Scalars<float, CPU, the_cpu_dev> Sc;
    typedef Vectors<float, CPU, the_cpu_dev> Vc;
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
    
    //Reading data
    bool readFromFile = true;
    OperatorsMats<Sc::value_type> mats(myIO, par.p_, 2*par.p_, readFromFile);
    
    // memory allocation
    Sur S(par, mats);  

    // initializing vesicle positions from text file
    char fname[400];
    sprintf(fname,"%s/precomputed/two_ellipse_12", 
        getenv("VES3D_DIR"));
    myIO.ReadData(fname,dLen,S.x_.begin());

    float T = 10;
    float dt = .001;
    CurvatureFlow<Sur> Cf;
    Cf(S, T, dt);
    
    myIO.WriteData("EvolveSurf.txt", S.x_.Size(), S.x_.begin());
}


