#include "Vectors.h"
#include "Surface.h"
#include "OperatorsMats.h"
#include "VesInteraction.h"
#include "InterfacialVelocity.h"

typedef double real;
extern const Device<CPU> the_cpu_device(0);

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Stokes Kernels Test:"
        <<"\n ==============================\n");

    typedef Scalars<real, CPU, the_cpu_device> Sca_t;
    typedef Vectors<real, CPU, the_cpu_device> Vec_t;

    int p = 12;
    int np = 2*p*(p+1);
    int nv = 100000;
    int rep = 1;
    ///@bug the number of targets effects the speedup dramatically
    int numtrg = 1;

    Vec_t pos(nv,p), den(nv,p), qw(1,p), pot1(nv,p), pot2(nv,p);
    
    pos.fillRand();
    den.fillRand();
    qw.fillRand();

    Logger::Tic();
    DirectStokesKernel(pos.getStride(), pos.getNumSubs(), 0,
        numtrg, qw.begin(), pos.begin(), 
        pos.begin(), den.begin(), pot1.begin());
    cout<<"  Time (Direct) : "<<Logger::Toc()<<endl;
    
    Logger::Tic();
    DirectStokesSSE(pos.getStride(), pos.getNumSubs(), 0,
        numtrg, qw.begin(), pos.begin(), 
        pos.begin(), den.begin(), pot2.begin());
    cout<<"  Time (SSE)    : "<<Logger::Toc()<<endl;
    
    axpy(-1,pot1,pot2, pot1);
    cout<<"\n  Error = "<<MaxAbs(pot1)<<endl;
}
