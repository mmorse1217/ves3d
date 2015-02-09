#include "PVFMMInterface.h"
#include "ves3d_common.h"
#include "Logger.h"
#include "EvolveSurface.h"
#include "ParallelLinSolver_Petsc.h"

typedef Device<CPU> Dev;
extern const Dev the_dev(0);

// Default callback for errors
Error_t cb_abort(const ErrorEvent &err)
{
    CERR_LOC("Aborting, received error "<<err,"",abort());
    return err.err_;
}

void run_sim(int argc, char **argv){
    SET_ERR_CALLBACK(&cb_abort);

    typedef EvolveSurface<real_t, Dev, the_dev> Evolve_t;
    typedef Evolve_t::Params_t Par_t;
    typedef Evolve_t::Arr_t Arr_t;
    typedef Evolve_t::Vec_t Vec_t;
    typedef Evolve_t::Interaction_t Inter_t;
    typedef ParallelLinSolverPetsc<real_t> PSol_t;

    // Setting the parameters
    Par_t sim_par;
    CHK(sim_par.parseInput(argc, argv));
    COUT(sim_par);

    //Initial vesicle positions
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);

    //reading the prototype form file
    DataIO myIO(FullPath(sim_par.save_file_name));

    char fname[300];
    std::string prec = (typeid(real_t) == typeid(float)) ? "float" : "double";
    sprintf(fname, sim_par.init_file_name.c_str(), sim_par.sh_order,prec.c_str());
    myIO.ReadData( FullPath(fname), x0, DataIO::ASCII, 0, x0.getSubLength());

    //reading centers file
    if (sim_par.cntrs_file_name.size()){
        INFO("Reading centers from file");
        Arr_t cntrs(DIM * sim_par.n_surfs);
        sprintf(fname, sim_par.cntrs_file_name.c_str(), sim_par.sh_order,prec.c_str());
        myIO.ReadData( FullPath(fname), cntrs, DataIO::ASCII, 0, cntrs.size());

        INFO("Populating the initial configuration using centers");
        Populate(x0, cntrs);
    };

    //Reading Operators From File
    bool readFromFile = true;
    Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

    // interaction handler
    Inter_t fmm_interaction(&PVFMMEval, &PVFMMDestroyContext<real_t>);
    // Inter_t fmm_interaction(&StokesAlltoAll); /* for debugging -- to be removed */

    // parallel solver
    PSol_t ksp(VES3D_COMM_WORLD);

    //Evolve surface class
    Evolve_t Es(sim_par, Mats, x0, &vInf, NULL, &fmm_interaction, NULL, &ksp);
    CHK( Es.Evolve() );
}

int main(int argc, char **argv)
{
    pvfmm::Profile::Enable(true);
    PROFILESTART();
    PetscInitialize(&argc, &argv, NULL, NULL);
    run_sim(argc, argv);
    PROFILEEND("",0);
    PRINTERRORLOG();
    PROFILEREPORT(SortTime);
    PetscFinalize();
    return 0;
}
