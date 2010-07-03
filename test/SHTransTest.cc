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
    typedef typename containers::Scalars<real,CPU,the_cpu_dev> Sca;
    typedef typename containers::Vectors<real,CPU,the_cpu_dev> Vec;
    
    int nVec(1), p(12), q(2*p);
    Parameters<real> params;
    params.sh_order = p;
    params.rep_up_freq = q;
    cout<<params<<endl;

    //reading mats

    DataIO<float, CPU> IO(the_cpu_dev);
    bool readFromFile(true);
    OperatorsMats<real, DataIO<float, CPU> > mats(IO, readFromFile, params);
    
    //Initializing vesicle positions from text file
    Vec x(nVec, p);
    char fname[200];
    sprintf(fname, "precomputed/biconcave_ra65_%d",p);
    IO.ReadData(fname, x.size(), x.begin());

    //Operators
    SHTrans<Sca> sht_p(p, mats.mats_p_, 1);
    SHTrans<Sca> sht_q(q, mats.mats_p_up_);
    
    //Work vectors
    int mx = (p > q)? p : q;
    Vec shc(nVec, mx), wrk(nVec, mx);
    
    //resampling
    Vec xr(nVec, q);
    Resample(x, sht_p, sht_q, shc, wrk, xr);
    IO.WriteData("ShtResampling.out", xr.size(), xr.begin(), ios::out);
    
    //Filtering
    IO.WriteData("ShtFiltering.out", x.size(), x.begin(), ios::out);
    sht_p.lowPassFilter(x, wrk, shc, x);
    IO.WriteData("ShtFiltering.out", x.size(), x.begin(), ios::app);


    return 0;
}
