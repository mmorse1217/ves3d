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

    bool readFromFile = true;
    Mats_t mats(the_cpu_device, myIO, readFromFile, sim_par);

    int np = x0.getStride();
    Sca_t rot_mat(p + 1, 1, make_pair(2 * p, np));
    Sca_t all_rot_mats(p + 1, 1, make_pair(np, np));
    Sca_t sp_harm_rot_mats(p+1, 1, make_pair(1, (p + 1) * (2 * p +  1) * (2 * p + 3)/3));

    fname = "precomputed/SpHarmRotMats_p12_float.bin";
    myIO.ReadData(fname, sp_harm_rot_mats);
    
    MovePole<Sca_t> move_pole(all_rot_mats, rot_mat, sp_harm_rot_mats);
    
    const Sca_t* inputs[] = {&x0};
    Sca_t* outputs[] = {&xr};
    move_pole.setOperands(inputs, 1, ViaSpHarm);
    move_pole(3, 0, outputs);

    myIO.Append(xr);    
}
