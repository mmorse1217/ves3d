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

typedef float real;
#define DT CPU

extern const Device<DT> the_device(0);
using namespace std;

template<typename Container>
void copyIdenticalTheFirst(Container &x)
{
    int cpysize = x.getSubN(1) - x.begin();
    cpysize *=sizeof(typename Container::value_type);

    for(int ii=1; ii<x.getNumSubs(); ++ii)
        x.getDevice().Memcpy(x.getSubN(ii), x.begin(), 
            cpysize, MemcpyDeviceToDevice);
}

int main(int argc, char** argv)
{
    typedef Scalars<real, DT, the_device> Sca_t;
    typedef Vectors<real, DT, the_device> Vec_t;
    typedef Surface<Sca_t,Vec_t> Sur_t;
    typedef VesInteraction<real> Interaction_t;
    typedef InterfacialVelocity<Sur_t, Interaction_t> IntVel_t;

    // Parameter
    int p(12), nvec(1);
    real tol(1e-5);
    Parameters<real> sim_par;
    sim_par.n_surfs = nvec;   
    sim_par.sh_order = p;
    sim_par.filter_freq = 8;
    sim_par.bending_modulus = 1e-2;
    sim_par.position_solver_tol = 1e-1 * tol;
    sim_par.tension_solver_tol = 1e-1 * tol;
    sim_par.ts = 1;
    sim_par.bg_flow_param = 1e-1;
    sim_par.singular_stokes = Direct;
    sim_par.upsample_interaction = false;//Make sure this is false

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
    OperatorsMats<Sca_t> Mats(readFromFile, sim_par);

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
    x.resize(nvec);       copyIdenticalTheFirst(x);
    Fb.resize(nvec);      copyIdenticalTheFirst(Fb);
    SFb.resize(nvec);     copyIdenticalTheFirst(SFb);
    vel.resize(nvec);     copyIdenticalTheFirst(vel);
    xnew.resize(nvec);    copyIdenticalTheFirst(xnew);
    tension.resize(nvec); copyIdenticalTheFirst(tension);
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
