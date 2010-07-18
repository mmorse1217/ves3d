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
    typedef OperatorsMats<Sca_t> Mats_t;
    int p(12);
    int nVec(1);
    
    Parameters<real> sim_par;
    char *fname = "MovePole.out";
    DataIO myIO(fname);
    remove(fname);
    
    Vec_t x0(nVec, p), xr(nVec, p);
    int fLen = x0.getStride();
  
    x0.fillRand();
    fname = "precomputed/dumbbell_cart12_single.txt";
    myIO.ReadData(fname, x0, 0, fLen * DIM);
    
    bool readFromFile = true;
    Mats_t mats(readFromFile, sim_par);

    SHTrans<Sca_t, SHTMats<real, Device<DT> > > sht(p, mats.mats_p_);
    MovePole<Sca_t, Mats_t> move_pole(mats);
        
    const Sca_t* inputs[] = {&x0};
    Sca_t* outputs[] = {&xr};
    
    //Correctness
    move_pole.setOperands(inputs, Direct);
    for(int ii=0; ii<x0.getGridDim().first; ++ii)
        for(int jj=0; jj<x0.getGridDim().second; ++jj)
        {    
            move_pole(ii, jj, outputs);
            myIO.Append(xr);    
        }

//     //Profile
//     nVec = 100;
//     x0.resize(nVec);
//     xr.resize(nVec);
//     int rep(2);
//     PROFILESTART();
//     move_pole.setOperands(inputs, Direct);
//     for(int kk=0;kk<rep; ++kk)
//         for(int ii=0; ii<x0.getGridDim().first; ++ii)
//             for(int jj=0; jj<x0.getGridDim().second; ++jj)
//             {    
//                 move_pole(ii, jj, outputs);
//                 //myIO.Append(xr);    
//             }
//     PROFILEEND("Direct rotation",0);

//     PROFILESTART();
//     move_pole.setOperands(inputs, ViaSpHarm);
//     for(int kk=0;kk<rep; ++kk)
//         for(int ii=0; ii<x0.getGridDim().first; ++ii)
//             for(int jj=0; jj<x0.getGridDim().second; ++jj)
//             {    
//                 move_pole(ii, jj, outputs);
//                 //myIO.Append(xr);    
//             }
//     PROFILEEND("SpHarm rotation",0);
//     PROFILEREPORT(SortTime);
}
