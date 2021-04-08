#include "ves3d_simulation.h"

typedef Device<CPU> Dev;
extern const Dev cpu(0);
typedef Simulation<Dev, cpu> Sim_t;
typedef Sim_t::Param_t Param_t;

int main(int argc, char **argv)
{
    SET_ERR_CALLBACK(&cb_abort);
    PROFILESTART();

    int pargc(0);
    char **pargv(NULL);
    // don't use petsc for commandline argument parsing
    VES3D_INITIALIZE(&pargc, &pargv, NULL, NULL);

    pvfmm::SetSigHandler();

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

    {
        // putting sim in block so that it is deconstructed before
        // petscfinalize
        Sim_t sim(argc,argv,&dict);
        CHK(sim.Run());
    }

    PROFILEEND("",0);
    PRINTERRORLOG();
    PROFILEREPORT(SortTime);
    VES3D_FINALIZE();
    return 0;
}
