#include <iostream>             
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "Parameters.h"
#include "MovePole.h"

#define DT CPU
typedef double real;
extern const Device<DT> the_cpu_device(0);


int main(int argc, char** argv)
{
    bool check_correct = (argc > 1) ? atoi(argv[1]) : true;
    bool profile = (argc > 2) ? atoi(argv[2]) : true;
    
    typedef Scalars<real, DT, the_cpu_device> Sca_t;
    typedef Vectors<real, DT, the_cpu_device> Vec_t;
    typedef OperatorsMats<Sca_t> Mats_t;
    int p(64);
    int nVec(1);
    
    Parameters<real> sim_par;
    sim_par.sh_order = p;
    char fname[300] = "MovePole.out";
    DataIO myIO(fname);
    remove(fname);
    
    Vec_t x0(nVec, p);
    int fLen = x0.getStride();
    string prec = (typeid(real) == typeid(float)) ? "float" : "double"; 
    sprintf(fname,"precomputed/dumbbell_%u_%s.txt",p,prec.c_str());
    myIO.ReadData(fname, x0, 0, fLen * DIM);
    
    bool readFromFile = true;
 //    sim_par.singular_stokes = Direct;
//     Mats_t mats_direct(readFromFile, sim_par);
    
    sim_par.singular_stokes = ViaSpHarm;
    Mats_t mats_spHarm(readFromFile, sim_par);

    SHTrans<Sca_t, SHTMats<real, Device<DT> > > sht(p, mats_spHarm.mats_p_);
    //MovePole<Sca_t, Mats_t> move_pole_direct(mats_direct);
    MovePole<Sca_t, Mats_t> move_pole_spHarm(mats_spHarm);
      
    const Sca_t* inputs[] = {&x0};
                
    //Correctness
    if ( check_correct )
    {
        Vec_t xr_direct(nVec, p), xr_spHarm(nVec, p);
        myIO.Append(x0);
        real err = 0;
        for(int ii=0; ii<x0.getGridDim().first; ++ii)
            for(int jj=0; jj<1; ++jj) //jj<x0.getGridDim().second; ++jj) 
            {
                // {
//                     move_pole_direct.setOperands(inputs, 1, Direct);
//                     Sca_t* output[] = {&xr_direct};
//                     move_pole_direct(ii, jj, output);
//                     myIO.Append(xr_direct);
//                 }
                
                
                {
                    move_pole_spHarm.setOperands(inputs, 1, ViaSpHarm);
                    Sca_t* output[] = {&xr_spHarm};
                    move_pole_spHarm(ii, jj, output);
                    myIO.Append(xr_spHarm);
                }
                axpy(-1,xr_direct, xr_spHarm, xr_spHarm);
                err = max(err,AlgebraicDot(xr_spHarm,xr_spHarm));
            }
        COUT("    The difference between the two methods : "<<err<<endl<<endl);
    }
    
//     //Profile
//     if ( profile )
//     {
//         nVec = 100;
//         x0.resize(nVec);
//         Vec_t xr(p,nVec);
//         int rep(2);
//         Sca_t* output[] = {&xr};

//         PROFILECLEAR();
//         PROFILESTART();
//         move_pole.setOperands(inputs, Direct);
//         for(int kk=0;kk<rep; ++kk)
//             for(int ii=0; ii<x0.getGridDim().first; ++ii)
//                 for(int jj=0; jj<x0.getGridDim().second; ++jj)
//                     move_pole(ii, jj, output);
    
//         PROFILEEND("Direct rotation ",0);
//         PROFILEREPORT(SortTime);
        
//         PROFILECLEAR();
//         PROFILESTART();
//         move_pole.setOperands(inputs, ViaSpHarm);
//         for(int kk=0;kk<rep; ++kk)
//             for(int ii=0; ii<x0.getGridDim().first; ++ii)
//                 for(int jj=0; jj<x0.getGridDim().second; ++jj)
//                     move_pole(ii, jj, output);
        
//         PROFILEEND("SpHarm rotation ",0);
//         PROFILEREPORT(SortTime);
//     }
}
