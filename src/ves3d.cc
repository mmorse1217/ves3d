#include "ves3d_simulation.h"

typedef Device<CPU> Dev;
extern const Dev cpu(0);
typedef Simulation<Dev, cpu> Sim_t;
typedef Sim_t::Param_t Param_t;

int main(int argc, char **argv)
{
    SET_ERR_CALLBACK(&cb_abort);
    PROFILESTART();

#ifdef HAS_PETSC
    PetscInitialize(&argc, &argv, NULL, NULL);
#endif

    DictString_t dict;
    int nproc(1), rank(0);
#ifdef HAS_MPI
    // Adding nproc and rank to template expansion dictionary
    MPI_Comm_size(VES3D_COMM_WORLD, &nproc);
    MPI_Comm_rank(VES3D_COMM_WORLD, &rank);
#endif
    std::stringstream snp, sr;
    snp<<nproc;
    sr<<std::setfill('0')<<std::setw(6)<<rank;
    dict["nprocs"] = snp.str();
    dict["rank"]   = sr.str();

    Param_t input_params;
    CHK(input_params.parseInput(argc, argv, &dict));

    {	// putting sim in block so that it is deconstructed before
	// petscfinalize
	Sim_t sim(input_params);
	CHK(sim.Run());
    }

    PROFILEEND("",0);
    PRINTERRORLOG();
    PROFILEREPORT(SortTime);

#ifdef HAS_PETSC
    PetscFinalize();
#endif

    return 0;
}
