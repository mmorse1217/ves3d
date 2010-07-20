#include "Vectors.h"
#include "Surface.h"
#include "OperatorsMats.h"
#include "VesInteraction.h"
#include "InterfacialVelocity.h"

typedef float real;

extern const Device<CPU> the_cpu_device(0);

template<typename Sca, typename Vec, typename Device>
void stokesTest(const Device &dev)
{
    typedef OperatorsMats<Sca> Mats_t;
    typedef Surface<Sca, Vec> Surf_t;
    typedef VesInteraction<real> Interaction_t;
    typedef InterfacialVelocity<Surf_t, Interaction_t> IntVel_t;

    int p(12);
    int nVec(1);
    
    Parameters<real> sim_par;
    DataIO myIO;
    
    Vec x(nVec, p);
    int fLen = x.getStride();
        
    char fname[400];
    sprintf(fname, "%s/precomputed/dumbbell_cart12_single.txt",
        getenv("VES3D_DIR"));
    myIO.ReadData(fname, x);
    
    bool readFromFile = true;
    Mats_t mats(readFromFile, sim_par);
        
    //Creating objects
    Surf_t S(x, mats);
    Interaction_t Interaction(NULL);
    IntVel_t intvel(S, Interaction, mats, sim_par);
    
    //Test solver
    BiCGStab<Vec, IntVel_t> solver;
    Vec xx(nVec,p), rhs(nVec, p), b(nVec,p);
    int mIter(100);
    real tol(1e-4);
    
    unsigned short tt = time(NULL);
    seed48(&tt);
    
    xx.fillRand();
    intvel(xx, rhs);

    axpy(static_cast<real>(0), xx, xx);

    solver(intvel, xx, rhs, mIter, tol);
    intvel(xx, b);
    axpy(static_cast<real>(-1),b,rhs,b);
    COUT( "\n\n                 relres = "<<tol
        <<"\n            True relres = "
        <<sqrt(AlgebraicDot(b,b)/AlgebraicDot(rhs,rhs))<<endl);
}

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Stokes Test:"
        <<"\n ==============================\n");
    sleep(1);

    typedef Scalars<real, CPU, the_cpu_device> Sca_t;
    typedef Vectors<real, CPU, the_cpu_device> Vec_t;
    
    //stokesTest<Sca_t, Vec_t, Device<CPU> >(the_cpu_device);
    
    int p = 12;
    int nv = 4000;
    int rep = 1;
    
    Vec_t pos(nv,p), den(nv,p), qw(1,p), pot1(nv,p), pot2(nv,p);
    
    pos.fillRand();
    den.fillRand();
    qw.fillRand();
    
    Logger::Tic();
    for(int ii=0;ii<rep; ++ii)
        DirectStokesKernel(pos.getStride(), pos.getNumSubs(), 0,
            pos.getStride(), qw.begin(), pos.begin(), 
            pos.begin(), den.begin(), pot1.begin());
    cout<<"  Direct :"<<Logger::Toc()<<endl;
    
    Logger::Tic();
    for(int ii=0;ii<rep; ++ii)
        DirectStokesSSE(pos.getStride(), pos.getNumSubs(), 0,
            pos.getStride(), qw.begin(), pos.begin(), 
            pos.begin(), den.begin(), pot2.begin());
    cout<<"  SSE    :"<<Logger::Toc()<<endl;
    
    axpy(-1,pot1,pot2, pot1);
    cout<<"\n  Error = "<<MaxAbs(pot1)<<endl;

    return 0;
}
