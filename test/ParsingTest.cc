#include "Parameters.h"
#include "Logger.h"

typedef double real;

int main(int argc, char **argv)
{
    VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
    Parameters<real>  sim_par(argc, argv);
    COUT(sim_par);
    VES3D_FINALIZE();
}
