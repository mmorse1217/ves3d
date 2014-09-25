#include "Parameters.h"
#include "Logger.h"

typedef double real;

int main(int argc, char **argv)
{
    Parameters<real>  sim_par(argc, argv);
    COUT(sim_par);
}
