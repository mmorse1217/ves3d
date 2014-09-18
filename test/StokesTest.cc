#include "Vectors.h"
#include "Surface.h"
#include "OperatorsMats.h"
#include "VesInteraction.h"
#include "InterfacialVelocity.h"
#include "HelperFuns.h"

typedef double real;
typedef Device<CPU> DevCPU;

extern const DevCPU the_cpu_dev(0);


int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Stokes Kernels Test:"
        <<"\n ==============================\n");

    typedef Scalars<real, DevCPU, the_cpu_dev> Sca_t;
    typedef Vectors<real, DevCPU, the_cpu_dev> Vec_t;

    int p = 12;
    int np = 2*p*(p+1);
    int nv = 100000;
    int rep = 1;
    ///@bug the number of targets effects the speedup dramatically
    int numtrg = 1;

    Vec_t pos(nv,p), den(nv,p), qw(1,p), pot1(nv,p), pot2(nv,p);

    fillRand(pos);
    fillRand(den);
    fillRand(qw);

    Logger::Tic();
    DirectStokesKernel(pos.getStride(), pos.getNumSubs(), 0,
        numtrg, qw.begin(), pos.begin(),
        pos.begin(), den.begin(), pot1.begin());
    COUT("  Time (Direct) : "<<Logger::Toc()<<std::endl);

    Logger::Tic();
    DirectStokesSSE(pos.getStride(), pos.getNumSubs(), 0,
        numtrg, qw.begin(), pos.begin(),
        pos.begin(), den.begin(), pot2.begin());
    COUT("  Time (SSE)    : "<<Logger::Toc()<<std::endl);

    axpy(-1,pot1,pot2, pot1);
    COUT("\n  Error = "<<MaxAbs(pot1)<<std::endl);
}
