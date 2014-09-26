#include <iostream>
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "Parameters.h"
#include "MovePole.h"

#define DT CPU
typedef double real;
typedef Device<DT> Dev;
extern const Dev the_device(0);

int main(int argc, char** argv)
{
    bool check_correct = (argc > 1) ? atoi(argv[1]) : true;
    bool profile = (argc > 2) ? atoi(argv[2]) : true;

    typedef Scalars<real, Dev, the_device> Sca_t;
    typedef Vectors<real, Dev, the_device> Vec_t;
    typedef typename Sca_t::array_type Arr;
    typedef OperatorsMats<Arr> Mats_t;

    int p(6);
    int nVec(1);

    Parameters<real> sim_par;
    sim_par.sh_order = p;
    sim_par.rep_up_freq = 2*p;

    char fname[300] = "MovePole.out";
    remove(fname);
    DataIO myIO(fname,DataIO::BIN);

    Vec_t x0(nVec, p);
    int fLen = x0.getStride();
    std::string prec = (typeid(real) == typeid(float)) ? "float" : "double";
    sprintf(fname,"precomputed/dumbbell_%u_%s.txt",p,prec.c_str());
    myIO.ReadData(fname, x0, DataIO::ASCII, 0, fLen * DIM);

    bool readFromFile = true;
    sim_par.singular_stokes = Direct;
    Mats_t mats_direct(readFromFile, sim_par);

    sim_par.singular_stokes = ViaSpHarm;
    Mats_t mats_spHarm(readFromFile, sim_par);

    MovePole<Sca_t, Mats_t> move_pole_direct(mats_direct);
    MovePole<Sca_t, Mats_t> move_pole_spHarm(mats_spHarm);

    const Sca_t* inputs[] = {&x0};

    // Correctness
    if ( check_correct )
    {
        Vec_t xr_direct(nVec, p), xr_spHarm(nVec, p);
        myIO.Append(x0);
        real err = 0;

        move_pole_direct.setOperands(inputs, 1, Direct /*or DirectEagerEval*/);
        move_pole_spHarm.setOperands(inputs, 1, ViaSpHarm);

        for(int ii=0; ii<x0.getGridDim().first; ++ii)
            for(int jj=0; jj<x0.getGridDim().second; ++jj)
            {
                {
                    Sca_t* output[] = {&xr_direct};
                    move_pole_direct(ii, jj, output);
                    myIO.Append(xr_direct);
                }

                if (!NO_SPARSE_MATVEC)
                {
                    Sca_t* output[] = {&xr_spHarm};
                    move_pole_spHarm(ii, jj, output);
                    myIO.Append(xr_spHarm);

                    axpy(-1,xr_direct, xr_spHarm, xr_spHarm);
                    err = std::max(err,AlgebraicDot(xr_spHarm,xr_spHarm));
                }
            }

        if (!NO_SPARSE_MATVEC)
            COUT("The difference between the two methods : "<<err);

        myIO.FlushBuffer<real>();
        COUT(alert<<" *** Use ../matlab/MovePoleTest.m to check the result ***"<<alert);
    }


    SHTrans<Sca_t, SHTMats<real, Device<DT> > > sht(p, mats_spHarm.mats_p_);
    //Profile
    if ( profile )
    {
        nVec = 100;
        x0.resize(nVec);
        Vec_t xr(p,nVec);
        int rep(2);
        Sca_t* output[] = {&xr};

        PROFILECLEAR();
        PROFILESTART();
        for(int kk=0;kk<rep; ++kk)
            for(int ii=0; ii<x0.getGridDim().first; ++ii)
                for(int jj=0; jj<x0.getGridDim().second; ++jj)
                    move_pole_direct(ii, jj, output);

        PROFILEEND("Direct rotation ",0);
        PROFILEREPORT(SortTime);

        if (!NO_SPARSE_MATVEC){
            PROFILECLEAR();
            PROFILESTART();
            for(int kk=0;kk<rep; ++kk)
                for(int ii=0; ii<x0.getGridDim().first; ++ii)
                    for(int jj=0; jj<x0.getGridDim().second; ++jj)
                        move_pole_spHarm(ii, jj, output);

            PROFILEEND("SpHarm rotation ",0);
            PROFILEREPORT(SortTime);
        }
    }
}
