#include "Device.h"
#include "Scalars.h"
#include "Parameters.h"
#include "OperatorsMats.h"
#include "SHTrans.h"
#include "Logger.h"

#define DT CPU
using namespace std;
typedef float real_t;
typedef Device<DT> Device_t;

extern Device_t the_device(0);

int main(int argc, char *argv[])
{
    typedef Scalars<real_t, DT, the_device> Sca_t;
    typedef OperatorsMats<Sca_t>  Mats_t;
    typedef SHTMats<real_t, Device_t> SHTMats_t;
    typedef SHTrans<Sca_t, SHTMats_t> SHT_t; 

    int ns(2000), p(16);
    Parameters<real_t> params;
    params.sh_order = p;
    params.n_surfs = ns;
    
    //Reading matrices form file
    bool readFromFile(true);
    Mats_t mats(readFromFile, params);
    
    //Test case
    Sca_t x(ns, p), shc(ns,p), wrk(ns,p);
    SHT_t sht(p, mats.mats_p_);
    
    //Calling methods
    sht.forward(x, wrk, shc);
    sht.backward(shc, wrk, x);
    sht.backward_du(shc, wrk, x);
    sht.backward_dv(shc, wrk, x);
    sht.backward_d2u(shc, wrk, x);
    sht.backward_d2v(shc, wrk, x);
    sht.backward_duv(shc, wrk, x);
    
    //Reporting results
    PROFILEREPORT(SortFunName);
    return 0;
}
