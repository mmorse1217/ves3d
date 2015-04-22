#include "ves3d_common.h"
#include "Logger.h"
#include "EvolveSurface.h"

#ifdef HAS_PETSC
#include "ParallelLinSolver_Petsc.h"
#endif

typedef Device<CPU> Dev;
extern const Dev the_dev(0);

// Default callback for errors
Error_t cb_abort(const ErrorEvent &err)
{
    CERR_LOC("Aborting, received error "<<err,"",abort());
    return err.err_;
}

void run_sim(int argc, char **argv){
    typedef EvolveSurface<real_t, Dev, the_dev> Evolve_t;
    typedef Evolve_t::Params_t Par_t;
    typedef Evolve_t::Arr_t Arr_t;
    typedef Evolve_t::Vec_t Vec_t;
    typedef Evolve_t::Interaction_t Inter_t;

    // Setting the parameters
    Par_t sim_par;
    CHK(sim_par.parseInput(argc, argv));
    omp_set_num_threads(sim_par.num_threads);
    COUT(sim_par);

    //Initial vesicle positions
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);

    //reading the prototype form file
    DataIO myIO(FullPath(sim_par.checkpoint_file_name));
    myIO.ReadData( FullPath(sim_par.init_file_name), x0, DataIO::ASCII, 0, x0.getSubLength());

    //reading centers file
    if (sim_par.cntrs_file_name.size()){
        INFO("Reading centers from file");
        Arr_t cntrs(DIM * sim_par.n_surfs);
        myIO.ReadData( FullPath(sim_par.cntrs_file_name), cntrs, DataIO::ASCII, 0, cntrs.size());

        INFO("Populating the initial configuration using centers");
        Populate(x0, cntrs);
    };

    //Reading Operators From File
    Evolve_t::Mats_t Mats(true /*readFromFile*/, sim_par);

    //Setting the background flow
    BgFlowBase<Vec_t> *vInf(NULL);
    CHK(BgFlowFactory(sim_par, &vInf));

#ifdef HAS_PETSC
    ParallelLinSolverPetsc<real_t> ksp(VES3D_COMM_WORLD);
    ParallelLinSolver<real_t> *pksp(&ksp);
#else
    ParallelLinSolver<real_t> *pksp(NULL);
#endif

    Inter_t direct_interaction(&StokesAlltoAll);

    //Evolve surface class
    Evolve_t Es(sim_par, Mats, x0, vInf, NULL, &direct_interaction, NULL, pksp);
    CHK( Es.Evolve() );

    delete vInf;
}

int main(int argc, char **argv)
{
#ifdef HAS_PETSC
    PetscInitialize(&argc, &argv, NULL, NULL);
#endif

    SET_ERR_CALLBACK(&cb_abort);
    PROFILESTART();
    run_sim(argc, argv);
    PROFILEEND("",0);
    PRINTERRORLOG();
    PROFILEREPORT(SortTime);

#ifdef HAS_PETSC
    PetscFinalize();
#endif
}
