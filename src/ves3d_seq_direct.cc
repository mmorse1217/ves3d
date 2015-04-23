#include "ves3d_simulation.h"

typedef Device<CPU> Dev;
extern const Dev cpu(0);
typedef Simulation<Dev, cpu> Sim_t;
typedef Sim_t::Param_t Param_t;

int main(int argc, char **argv)
{
#ifdef HAS_PETSC
    PetscInitialize(&argc, &argv, NULL, NULL);
#endif

    SET_ERR_CALLBACK(&cb_abort);
    PROFILESTART();

    Param_t input_params;
    CHK(input_params.parseInput(argc, argv));
    Sim_t sim(input_params);
    CHK(sim.Run());

    PROFILEEND("",0);
    PRINTERRORLOG();
    PROFILEREPORT(SortTime);

#ifdef HAS_PETSC
    PetscFinalize();
#endif
}
