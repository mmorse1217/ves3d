#include <iostream>             
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "Parameters.h"
#include "MovePole.h"

#define DT CPU
typedef float real;
extern const Device<DT> the_cpu_device(0);


int main(int , char** )
{
    typedef Scalars<real, DT, the_cpu_device> Sca_t;
    typedef Vectors<real, DT, the_cpu_device> Vec_t;
    typedef OperatorsMats<real, Device<DT> > Mats_t;
    int const p(12);
    int const nVec(1);
    
    Parameters<real> sim_par;
    char *fname = "MovePole.out";
    DataIO myIO(fname);
    remove(fname);
    
    Vec_t x0(nVec, p), xr(nVec, p);
    int fLen = x0.getStride();
        
    fname = "precomputed/dumbbell_cart12_single.txt";
    myIO.ReadData(fname, fLen * x0.getTheDim(), x0.begin());
    Sca_t sp_harm_rot_mats(p+1, 1, make_pair(1, p * (4 * p * p -  1)/3 + 4 * p * p));  
    fname = "precomputed/SpHarmRotMats_p12_float.bin";
    myIO.ReadData(fname, sp_harm_rot_mats);

    bool readFromFile = true;
    Mats_t mats(the_cpu_device, myIO, readFromFile, sim_par);

    SHTrans<Sca_t, SHTMats<real, Device<DT> > > sht(p, mats.mats_p_);
    MovePole<Sca_t, Mats_t> move_pole(mats, sp_harm_rot_mats);
    
    Vec_t shc(nVec,p), wrk(nVec,p);
    
    const Sca_t* inputs[] = {&x0};
    Sca_t* outputs[] = {&xr};
    
    move_pole.setOperands(inputs, 1, ViaSpHarm);
    move_pole(3, 0, outputs);
    myIO.Append(xr);    
}
