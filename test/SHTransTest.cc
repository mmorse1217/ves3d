///This test only performs the filtering. Use matlab/SHTransTest.m to
///see the result in Matlab

#include "OperatorsMats.h"
#include "Device.h"
#include "DataIO.h"
#include "Parameters.h"
#include "HelperFuns.h"
#include "Scalars.h"
#include "Vectors.h"

typedef float real;

typedef Device<CPU> DCPU;
extern const DCPU the_cpu_dev(0);

int main(int argc, char *argv[])
{
    typedef Scalars<real,DCPU,the_cpu_dev> Sca_t;
    typedef typename Sca_t::array_type Arr_t;
    typedef Vectors<real,DCPU,the_cpu_dev> Vec_t;
    typedef OperatorsMats<Arr_t> OMats_t;
    typedef SHTMats<real,DCPU> SMats_t;
    typedef SHTrans<Sca_t,SMats_t> Sh_t;

    int nVec(1), p(6), q(2*p);
    Parameters<real> params;
    params.sh_order = p;
    params.rep_up_freq = q;
    COUT(params);

    //reading mats2
    DataIO IO;

    //Initializing vesicle positions from text file
    Vec_t x(nVec, p);
    char fname[200];
    sprintf(fname, "precomputed/dumbbell_%d_double.txt",p);
    IO.ReadData(fname, x);
    IO.WriteData("Original.out", x, std::ios::out);

    //Operators
    bool readFromFile(true);
    OMats_t M(readFromFile, params);

    //Work vectors
    int mx = (p > q)? p : q;
    Vec_t shc(nVec, mx), wrk(nVec, mx);

    //resampling
    Vec_t xr(nVec, q);
    Sh_t SH_p(params.sh_order,M.mats_p_);
    Sh_t SH_q(params.rep_up_freq,M.mats_p_up_);

    Resample(x, SH_p, SH_q, shc, wrk, xr);
    IO.WriteData("ShtResampling.out", xr, std::ios::out);

    //Filtering
    IO.WriteData("ShtFiltering.out", x, std::ios::out);
    SH_p.lowPassFilter(x, wrk, shc, x);
    IO.WriteData("ShtFiltering.out", x, std::ios::app);

    COUT(alert<<"\n  *** Run ../matlab/SHTransTest.m to see the results ***"
        <<alert);
    return 0;
}
