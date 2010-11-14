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
    int nSnapShots(100);
    int n_surfs(2), sh_order(12), n_eval(50 * 20 * 20);
    real bending_modulus(.01);
    Par_t sim_par;
    
    sim_par.sh_order = sh_order;
    sim_par.n_surfs = n_surfs;   
    sim_par.bending_modulus = bending_modulus;
    sim_par.bg_flow_param = 0.1;

    COUT(sim_par<<endl);

    //Building containers
    Sca_t saved_data( (DIM+1)* n_surfs * nSnapShots, sh_order);
    Vec_t x0(n_surfs, sh_order);
    Sca_t tension(n_surfs, sh_order);
    Vec_t x_eval(1, 1, make_pair(n_eval,1));
    Vec_t vel(1, 1, make_pair(n_eval,1));

    //Reading the prototype form file  
    DataIO myIO;
    char fname[] = "../matlab/StreamlineData.out";
    myIO.ReadData(fname, saved_data, 0, saved_data.size());

    //The evaluation points
    sprintf(fname,"../matlab/EvaluationPoints.out");
    myIO.ReadData(fname, x_eval, 0, x_eval.size());

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

    //Reading Operators From File
    Mats_t mats(true, sim_par);

    //Getting the velocity field
    Eval_t EV(&StokesAlltoAll, vInf, mats, bending_modulus);

    //Evaluating snapshots
    for ( int ii=0; ii<nSnapShots; ++ii )
    {       
        //Copying position and tension to their containers
        Sca_t::getDevice().Memcpy(x0.begin(), 
            saved_data.getSubN((DIM+1) * n_surfs * ii),
            x0.size() * sizeof(real), MemcpyDeviceToDevice);

        Sca_t::getDevice().Memcpy(tension.begin(), 
            saved_data.getSubN((DIM+1) * n_surfs * ii + DIM * n_surfs),
            tension.size() * sizeof(real), MemcpyDeviceToDevice);

        QC( EV(x0, tension, x_eval, vel) );
    
        //Writing to file
        myIO.WriteData("../matlab/Surfaces.out", x0, ios::app);
        myIO.WriteData("VelocityField.out", vel, ios::app);
    }

    return 0;
}
