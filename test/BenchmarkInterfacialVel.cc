#ifdef NDEBUG
#warning "This is a test file and would not compile in the non-debug mode. Remove the -DNDEBUG option and recompile."
#endif

#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
#include "OperatorsMats.h"
#include "InterfacialVelocity.h"
#include "BgFlow.h"
#include "CPUKernels.h"

#define DT CPU
typedef float real;
typedef Device<DT> Dev;

extern const Dev the_device(0);

using namespace std;

template<typename Container>
void populateByFirst(Container &x)
{
    int cpysize = x.getSubN_end(0) - x.getSubN_begin(0);
    cpysize *=sizeof(typename Container::value_type);

    for(int ii=1; ii<x.getNumSubs(); ++ii)
        x.getDevice().Memcpy(x.getSubN_begin(ii), x.begin(),
            cpysize, x.getDevice().MemcpyDeviceToDevice);
}

int main(int argc, char** argv)
{
    typedef Scalars<real, Dev, the_device> Sca_t;
    typedef Vectors<real, Dev, the_device> Vec_t;
    typedef typename Sca_t::array_type Arr_t;
    typedef OperatorsMats<Arr_t> Mats_t;
    typedef Surface<Sca_t,Vec_t> Sur_t;
    typedef VesInteraction<real> Interaction_t;
    typedef InterfacialVelocity<Sur_t, Interaction_t> IntVel_t;

    // Parameter
    int p(12), nvec(1);
    real tol(4e-5);
    Parameters<real> sim_par;
    sim_par.n_surfs              = nvec;
    sim_par.sh_order             = p;
    sim_par.filter_freq          = 8;
    sim_par.bending_modulus      = 1e-2;
    sim_par.position_solver_tol  = 1e-1 * tol;
    sim_par.tension_solver_tol   = 1e-1 *  tol;
    sim_par.ts                   = 1;
    sim_par.bg_flow_param        = 1e-1;
    sim_par.singular_stokes      = Direct;
    sim_par.upsample_interaction = false; //Make sure this is false
    sim_par.rep_up_freq          = 12;

    DataIO myIO;
    Vec_t x(nvec, p), Fb(nvec, p), SFb(nvec, p), vel(nvec, p), xnew(nvec, p);
    Sca_t tension(nvec, p);

    //reading the prototype form file
    myIO.ReadData("precomputed/Bmark_SingleVes_Exp_p12_float_x0.txt",x);
    myIO.ReadData("precomputed/Bmark_SingleVes_Exp_p12_float_Fb.txt",Fb);
    myIO.ReadData("precomputed/Bmark_SingleVes_Exp_p12_float_SFb.txt",SFb);
    myIO.ReadData("precomputed/Bmark_SingleVes_Exp_p12_float_vel.txt",vel);
    myIO.ReadData("precomputed/Bmark_SingleVes_Exp_p12_float_tension.txt",tension);
    myIO.ReadData("precomputed/Bmark_SingleVes_Exp_p12_float_xnew.txt",xnew);

    //Reading Operators From File
    bool readFromFile = true;
    Mats_t Mats(readFromFile, sim_par);

    //Making The Surface, And Time Stepper
    Sur_t S(x, Mats);

    {
        //Setting the background flow
        ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

        Interaction_t interaction(NULL);
        IntVel_t F(S, interaction, Mats, sim_par, vInf);

        // b-marking
        F.benchmarkExplicit(Fb, SFb, vel, tension, xnew, tol);
    }

    // Multiple vesicles no interaction
    nvec = 50;
    x.resize(nvec);       populateByFirst(x);
    Fb.resize(nvec);      populateByFirst(Fb);
    SFb.resize(nvec);     populateByFirst(SFb);
    vel.resize(nvec);     populateByFirst(vel);
    xnew.resize(nvec);    populateByFirst(xnew);
    tension.resize(nvec); populateByFirst(tension);
    S.setPosition(x);

    {
        //Setting the background flow
        ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

        Interaction_t Interaction(NULL);
        IntVel_t F(S, Interaction, Mats, sim_par,vInf);
        // b-marking
        F.benchmarkExplicit(Fb, SFb, vel, tension, xnew, tol);
    }

    // Multiple vesicles with interaction
    nvec = 2;
    x.resize(nvec);
    Fb.resize(nvec);
    SFb.resize(nvec);
    vel.resize(nvec);
    xnew.resize(nvec);
    tension.resize(nvec);

    myIO.ReadData("precomputed/Bmark_TwoVes_Exp_p12_float_x0.txt",x);
    myIO.ReadData("precomputed/Bmark_TwoVes_Exp_p12_float_Fb.txt",Fb);
    myIO.ReadData("precomputed/Bmark_TwoVes_Exp_p12_float_SFb.txt",SFb);
    myIO.ReadData("precomputed/Bmark_TwoVes_Exp_p12_float_vel.txt",vel);
    myIO.ReadData("precomputed/Bmark_TwoVes_Exp_p12_float_tension.txt",tension);
    myIO.ReadData("precomputed/Bmark_TwoVes_Exp_p12_float_xnew.txt",xnew);

    S.setPosition(x);
    {
        //Setting the background flow
        ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

        Interaction_t Interaction(StokesAlltoAll);
        IntVel_t F(S, Interaction, Mats, sim_par,vInf);
        // b-marking
        F.benchmarkExplicit(Fb, SFb, vel, tension, xnew, tol);
    }

    // Multiple vesicles implicit
    myIO.ReadData("precomputed/Bmark_TwoVes_Imp_p12_float_x0.txt",x);
    myIO.ReadData("precomputed/Bmark_TwoVes_Imp_p12_float_matvec.txt",vel);
    myIO.ReadData("precomputed/Bmark_TwoVes_Imp_p12_float_tension.txt",tension);
    myIO.ReadData("precomputed/Bmark_TwoVes_Imp_p12_float_xnew.txt",xnew);

    S.setPosition(x);
    {
        //Setting the background flow
        ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

        Interaction_t Interaction(StokesAlltoAll);
        IntVel_t F(S, Interaction, Mats, sim_par,vInf);
        // b-marking
        F.benchmarkImplicit(tension, vel, xnew, tol);
    }
}
