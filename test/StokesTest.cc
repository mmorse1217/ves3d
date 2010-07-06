#include "Vectors.h"
#include "Surface.h"
#include "OperatorsMats.h"
#include "VesInteraction.h"
#include "InterfacialVelocity.h"

typedef double real;

extern const Device<CPU> the_cpu_device(0);

template<typename Sca, typename Vec, typename Device>
void stokesTest(const Device &dev)
{
    typedef OperatorsMats<real, Device> Mats_t;
    typedef Surface<Sca, Vec> Surf_t;
    typedef VesInteraction<real> Interaction_t;
    typedef InterfacialVelocity<Surf_t, Interaction_t> IntVel_t;

    int p(12);
    int nVec(1);
    
    Parameters<real> sim_par;
    DataIO<real, Device> myIO(dev);
    
    Vec x(nVec, p);
    int fLen = x.getStride();
        
    char fname[400];
    sprintf(fname, "%s/precomputed/dumbbell_cart12_single.txt",
        getenv("VES3D_DIR"));
    myIO.ReadData(fname, fLen * x.getTheDim(), x.begin());
    
    bool readFromFile = true;
    Mats_t mats(dev, myIO, readFromFile, sim_par);
        
    //Creating objects
    Surf_t S(x, mats);
    Interaction_t Interaction(NULL);
    IntVel_t intvel(S, Interaction, mats, sim_par);
    
    //Test solver
    typedef Vec Cont;
    BiCGStab<Cont, IntVel_t> solver;
    Cont xx(nVec,p), rhs(nVec, p), b(nVec,p);
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

    typedef Scalars<real, CPU, the_cpu_device> ScaCPU_t;
    typedef Vectors<real, CPU, the_cpu_device> VecCPU_t;
    
    stokesTest<ScaCPU_t, VecCPU_t, Device<CPU> >(the_cpu_device);
    
    return 0;
}
