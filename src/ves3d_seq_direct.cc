#include "ves3d_common.h"
#include "Logger.h"
#include "EvolveSurface.h"

typedef double real;
typedef Device<CPU> Dev;
extern const Dev the_dev(0);

// Default callback for errors
Error_t cb_abort(const Error_t &err)
{
    CERR_LOC("Aborting, received error "<<err,"",abort());
    return err;
}

int main(int argc, char **argv)
{
    SET_ERR_CALLBACK(&cb_abort);

    typedef EvolveSurface<real, Dev, the_dev> Evolve_t;
    typedef Evolve_t::Params_t Par_t;
    typedef Evolve_t::Vec_t Vec_t;


    // Setting the parameters
    Par_t sim_par;
    CHK(sim_par.parseInput(argc, argv));
    COUT(sim_par);

    // //Initial vesicle positions
    // Vec_t x0(sim_par.n_surfs, sim_par.sh_order);

    // //reading the prototype form file
    // DataIO myIO(sim_par.save_file_name);
    // char fname[300];
    // std::string prec = (typeid(real) == typeid(float)) ? "float" : "double";
    // sprintf(fname,"%s/precomputed/dumbbell_%u_%s.txt", VES3D_PATH,
    //     sim_par.sh_order,prec.c_str());
    // myIO.ReadData(fname, x0, DataIO::ASCII, 0, x0.getSubLength());

    // //Reading Operators From File
    // bool readFromFile = true;
    // Evolve_t::Mats_t Mats(readFromFile, sim_par);

    // //Setting the background flow
    // ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

    // Evolve_t Es(sim_par, Mats, x0, &vInf);

    // CLEARERRORHIST();
    // PROFILESTART();
    // CHK( Es.Evolve() );
    // PROFILEEND("",0);

    // PRINTERRORLOG();
    // PROFILEREPORT(SortFlop);
    // myIO.Append(Es.S_->getPosition());
    // myIO.FlushBuffer<real>();
}
