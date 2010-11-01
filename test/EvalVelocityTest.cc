#include "DataIO.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "EvalVelocity.h"
#include "BgFlow.h"
#include "CPUKernels.h"

#include "Surface.h"

#define DT CPU
extern const Device<DT> the_device(0);
typedef float real;

int main(int argc, char** argv)
{
    typedef Parameters<real> Par_t;
    typedef Scalars<real, DT, the_device> Sca_t;
    typedef Vectors<real, DT, the_device> Vec_t;
    typedef OperatorsMats<Sca_t> Mats_t;
    typedef void(*StokesEval_t)(const real*, const real*, 
        size_t, real*, void*);    
    typedef EvalVelocity<Sca_t, Vec_t, StokesEval_t> Eval_t;

    //Parameters
    int n_surfs(1), sh_order(12), n_eval(10 * 10 * 10);
    real bending_modulus(0);
    Par_t sim_par;
    
    sim_par.sh_order = sh_order;
    sim_par.n_surfs = n_surfs;   
    sim_par.bending_modulus = bending_modulus;
    sim_par.bg_flow_param = 100;

    COUT(sim_par<<endl);

    //Building containers
    Vec_t x0(n_surfs, sh_order);
    Sca_t tension(n_surfs, sh_order);
    Vec_t x_eval(1, 1, make_pair(n_eval,1));
    Vec_t vel(1, 1, make_pair(n_eval,1));

    //Reading the prototype form file  
    DataIO myIO;
    char fname[300];
    sprintf(fname,"precomputed/dumbbell_%u_float.txt",sh_order);
    myIO.ReadData(fname, x0, 0, x0.size());

    sprintf(fname,"../matlab/EvaluationPoints.out");
    myIO.ReadData(fname, x_eval, 0, x_eval.size());

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

    //Reading Operators From File
    Mats_t mats(true, sim_par);

    //Getting the velocity field
    axpy(0,tension,tension);

    Eval_t EV(&StokesAlltoAll, vInf, mats, bending_modulus);
    QC( EV(x0, tension, x_eval, vel) );
    
    //Writing to file
    myIO.WriteData("../matlab/Surfaces.out", x0);
    myIO.WriteData("../matlab/EvaluationPoints.out", x_eval);
    myIO.WriteData("VelocityField.out", vel);

    return 0;
}
