///This test only performs the filtering. Use matlab/SHTransTest.m to
///see the result in Matlab

#include "OperatorsMats.h"
#include "Device.h"
#include "DataIO.h"
#include "Parameters.h"
#include "HelperFuns.h"
#include "Scalars.h"
#include "Vectors.h"

typedef double real;

typedef Device<CPU> DCPU;
extern const DCPU the_cpu_dev(0);

bool test_resample(){

    typedef Scalars<real,DCPU,the_cpu_dev> Sca_t;
    typedef typename Sca_t::array_type Arr_t;
    typedef Vectors<real,DCPU,the_cpu_dev> Vec_t;
    typedef OperatorsMats<Arr_t> OMats_t;
    typedef SHTMats<real,DCPU> SMats_t;
    typedef SHTrans<Sca_t,SMats_t> Sh_t;

    int nVec(1), p(6), q(2*p);
    Parameters<real> params;
    params.sh_order = p;
    params.upsample_freq = q;
    COUT(params);
    DataIO IO;

    //Initializing vesicle positions from text file
    Vec_t x(nVec, p);
    char fname[200];
    sprintf(fname, "precomputed/dumbbell_%d_double.txt",p);
    IO.ReadData(FullPath(fname), x, DataIO::ASCII);
    IO.WriteData("Original.out", x, DataIO::ASCII, std::ios::out);

    //Operators
    bool readFromFile(true);
    OMats_t M(readFromFile, params);

    //Work vectors
    int mx = (p > q)? p : q;
    Vec_t shc(nVec, mx), wrk(nVec, mx);

    //resampling
    Vec_t xr(nVec, q);
    Sh_t SH_p(params.sh_order,M.mats_p_);
    Sh_t SH_q(params.upsample_freq,M.mats_p_up_);

    Resample(x, SH_p, SH_q, shc, wrk, xr);
    IO.WriteData("ShtResampling.out", xr, DataIO::ASCII, std::ios::out);

    //Filtering
    IO.WriteData("ShtFiltering.out", x, DataIO::ASCII, std::ios::out);
    SH_p.lowPassFilter(x, wrk, shc, x);
    IO.WriteData("ShtFiltering.out", x, DataIO::ASCII, std::ios::app);

    COUT(alert<<"\n  *** Run ../matlab/SHTransTest.m to see the results ***"
        <<alert);

    return true;
}

bool test_inverse(){
    typedef Scalars<real,DCPU,the_cpu_dev> Sca_t;
    typedef typename Sca_t::array_type Arr_t;
    typedef OperatorsMats<Arr_t> OMats_t;
    typedef SHTMats<real,DCPU> SMats_t;
    typedef SHTrans<Sca_t,SMats_t> Sh_t;

    int p(6), q(2*p), n(1);
    Parameters<real> params;
    params.sh_order	 = p;
    params.upsample_freq = q;
    DataIO IO;

    // Operators
    OMats_t M(true /* readFromFile */, params);
    Sh_t sht(params.sh_order,M.mats_p_);
    Sca_t xhat(n,p), x(n,p), shc(n,p), wrk(n,p);
    size_t len((p+1)*(p+1)-1);

    fillRand(xhat);
    xhat.getDevice().Memset(xhat.begin()+len,0, (xhat.size()-len)*sizeof(real));

    sht.backward(xhat,wrk,x);
    sht.forward(x    ,wrk,shc);
    shc.getDevice().Memset(shc.begin()+len,0, (shc.size()-len)*sizeof(real));

    axpy(-1.0, shc,xhat,shc);
    real merr = MaxAbs(shc);
    ASSERT(merr<1e-14,"forward and backward shoud be inverse of each other, error="<<merr);
    return true;
}

int main(int argc, char *argv[])
{
    VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
    ASSERT(test_resample(),"resample test failed");
    ASSERT(test_inverse(),"inverse test failed");
    return 0;
    VES3D_FINALIZE();
}
