#include "EvolveSurface.h"
#include "CPUKernels.h"

#define DT CPU 
typedef double real;
extern const Device<DT> the_device(0);

int main(int argc, char **argv)
{
    typedef Parameters<real> Par_t;
    typedef EvolveSurface<real, DT, the_device> Evolve_t;
    typedef Evolve_t::Sca_t Sca_t;
    typedef Evolve_t::Vec_t Vec_t;

    // Setting the parameters
    Par_t sim_par(argc, argv);
    remove(sim_par.save_file_name.c_str());

    //Initial vesicle positions 
    Vec_t x0(sim_par.n_surfs, sim_par.sh_order);
    
    //reading the prototype form file  
    DataIO myIO(sim_par.save_file_name);
    myIO.ReadData(sim_par.init_file_name, x0, 0, x0.getSubLength());

    //Reading the Centers and populating
    if( !sim_par.cntrs_file_name.empty() )
    {
        Array<real, DT, the_device> cntrs(DIM * sim_par.n_surfs);
        myIO.ReadData(sim_par.cntrs_file_name, cntrs, 0, cntrs.size());
        Populate(x0, cntrs);
    }
    
    //Reading Operators From File
    bool readFromFile = true;
    Evolve_t::Mats_t Mats(readFromFile, sim_par);

    //Setting the background flow
    ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);
                
    Evolve_t Es(sim_par, Mats, x0, &vInf, &StokesAlltoAll);
        
    CLEARERRORHIST();
    PROFILESTART(); 
    QC( Es.Evolve() );    
    PROFILEEND("",0);
    
    PRINTERRORLOG();
    PROFILEREPORT(SortFlop);
}
