#include <iostream>             
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "Parameters.h"
#include "MovePole.h"

#define DT CPU
typedef float real;
extern const Device<DT> the_device(0);

int main(int , char** )
{
    typedef Scalars<real, DT, the_device> Sca_t;
    typedef Vectors<real, DT, the_device> Vec_t;
    typedef OperatorsMats<Sca_t> Mats_t;

    int nVec(1024);
    int P[] = {8, 12, 16, 24};
    int pLen = 4;
    for(int ii=0;ii<pLen; ++ii)
    {
        Parameters<real> sim_par;
        sim_par.sh_order = P[ii];
        sim_par.singular_stokes = Direct; //The files are read
                                          //according to this: Direct,
                                          //ViaSpHarm, DirectEagerEval
                                         
        Vec_t x0(nVec, P[ii]), xr(nVec, P[ii]);
        int fLen = x0.getStride();
        bool readFromFile = true;
        Mats_t mats(readFromFile, sim_par);
        
        SHTrans<Sca_t, SHTMats<real, Device<DT> > > sht(P[ii], mats.mats_p_);
        MovePole<Sca_t, Mats_t> move_pole(mats);
        
        const Sca_t* inputs[] = {&x0};
        Sca_t* outputs[] = {&xr};
    
        //Profile
        PROFILECLEAR();
        PROFILESTART();
        move_pole.setOperands(inputs, 1, sim_par.singular_stokes);
        for(int ii=0; ii<x0.getGridDim().first; ++ii)
            for(int jj=0; jj<x0.getGridDim().second; ++jj)
                move_pole(ii, jj, outputs);

        PROFILEEND("Direct_",0);
        PROFILEREPORT(SortTime);
    }
}
