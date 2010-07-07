#include "Vectors.h"
#include "Parameters.h"
#include "Surface.h"
#include "VesInteraction.h"
#include "OperatorsMats.h"
#include "InterfacialVelocity.h"

typedef float real;
#define DT CPU

extern const Device<DT> the_device(0);

int main(int argc, char** argv)
{
    typedef Scalars<real, DT, the_device> Sca_t;
    typedef Vectors<real, DT, the_device> Vec_t;
    typedef Surface<Sca_t,Vec_t> Sur_t;
    typedef VesInteraction<real> Interaction_t;
    typedef InterfacialVelocity<Sur_t, Interaction_t> IntVel_t;

    Parameters<real> sim_par; cout<<sim_par<<endl;
    DataIO<real, Device<DT> > myIO(the_device);

    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",x0);
    
    // The Interaction Class
    Interaction_t Interaction(NULL);

    //Reading Operators From File
    bool readFromFile = true;
    OperatorsMats<real, Device<DT> > Mats(the_device, myIO, 
        readFromFile, sim_par);

    //Making The Surface, And Time Stepper
    Sur_t S(x0, Mats);
    IntVel_t F(S, Interaction, Mats, sim_par);

    real tol = 1e-3;
    Sca_t tmp;
    
    //benchmarkExplicit(Vec_t &Fb, Vec_t &SFb, Sca_t &tension, Vec_t &vel, Vec_t &xnew, value_type tol) const
    F.benchmarkExplicit(x0, x0, tmp, x0, x0, tol);
}
