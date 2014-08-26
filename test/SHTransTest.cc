///This test only performs the filtering. Use matlab/SHTransTest.m to
///see the result in Matlab

#include "Device.h"
#include "OperatorsMats.h"
#include "DataIO.h"
#include "Parameters.h"
#include "Surface.h"
#include "Scalars.h"
#include "Vectors.h"

using namespace std;
typedef float real;

extern Device<CPU> the_cpu_dev(0);

int main(int argc, char *argv[])
{
    typedef Scalars<real,CPU,the_cpu_dev> Sca_t;
    typedef Vectors<real,CPU,the_cpu_dev> Vec_t;
    typedef OperatorsMats<Sca_t> Mats_t;
    typedef SHTrans<Sca_t,Mats_t> Sh_t;

    int nVec(1), p(6), q(2*p);
    Parameters<real> params;
    params.sh_order = p;
    params.rep_up_freq = q;
    cout<<params<<endl;

    //reading mats
    DataIO IO;

    //Initializing vesicle positions from text file
    Vec_t x(nVec, p);
    char fname[200];
    sprintf(fname, "precomputed/biconcave_ra65_%d",p);
    IO.ReadData(fname, x);

    //Operators
    bool readFromFile(true);
    Mats_t M(readFromFile, params);

    //Work vectors
    int mx = (p > q)? p : q;
    Vec_t shc(nVec, mx), wrk(nVec, mx);

    //resampling
    Vec_t xr(nVec, q);
    Sh_t SH_p(params.sh_order,M);
    Sh_t SH_q(params.rep_up_freq,M);
    assert(false); //outdated test
    //Resample(x, SH_p, SH_q, shc, wrk, xr);
    IO.WriteData("ShtResampling.out", xr, ios::out);

    // //Filtering

    // IO.WriteData("ShtFiltering.out", x, ios::out);
    // SH.lowPassFilter(x, wrk, shc, x);
    // IO.WriteData("ShtFiltering.out", x, ios::app);

    return 0;
}
