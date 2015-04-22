#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "DataIO.h"
#include "Parameters.h"
#include "Surface.h"
#include "Repartition.h"

typedef float real;
#define DT CPU

typedef Device<DT> Dev;
extern const Dev the_device(0);

#ifndef Doxygen_skip

// assigning more surfaces
void MyRepart(size_t nv, size_t stride, const real* x,
    const real* tension, size_t* nvr,
    real** xr, real** tensionr, void**)
{
    COUTDEBUG(emph<<"master repartitioning"<<emph);
    *nvr = 2*nv;
    *xr = new real[*nvr * stride * DIM];
    *tensionr = new real[*nvr * stride];

    int length = *nvr * DIM * stride;
    for(int ii=0; ii<length;++ii)
        *(*xr +ii) = ii;

    length = *nvr * stride;
    for(int ii=0; ii<length;++ii)
        *(*tensionr + ii) = ii;
}
#endif

int main(int argc, char** argv)
{
    //Checking that changing the number of surfaces at run-time works
    typedef Scalars<real, Dev, the_device> Sca_t;
    typedef Vectors<real, Dev, the_device> Vec_t;
    typedef typename Sca_t::array_type Arr_t;
    typedef OperatorsMats<Arr_t> Mats_t;

    {
        int const p(6);
        int const nv1(1), nv2(4), nv3(4000);

        //Reading operators from file
        Parameters<real> sim_par;
        sim_par.sh_order = p;
        sim_par.rep_up_freq = p;
        bool readFromFile = true;
        Mats_t mats(readFromFile, sim_par);

        //Creating objects
        Vec_t x1(nv1, p); fillRand(x1);
        Vec_t x2(nv2, p); fillRand(x2);
        Vec_t x3(nv3, p); fillRand(x3);

        Surface<Sca_t, Vec_t> S(mats,&x1);
        COUT("Number of surfaces: "<<S.getNumberOfSurfaces());
        Sca_t H1(nv1, p);
        axpy(1, S.getMeanCurv(), H1);
        ASSERT(S.getNumberOfSurfaces()==nv1,"Expeced number of surfaces");

        //Increasing the number of vesicles
        S.getPositionModifiable().resize(nv2);
        ASSERT(S.getNumberOfSurfaces()==nv2,"Expeced number of surfaces");
        axpy(1, x2, S.getPositionModifiable());
        Sca_t H2(nv2, p);
        axpy(1, S.getMeanCurv(), H2);
        COUT("Number of surfaces: "<<S.getNumberOfSurfaces());

        //Increasing the number of vesicles
        S.getPositionModifiable().resize(nv3);
        ASSERT(S.getNumberOfSurfaces()==nv3,"Expeced number of surfaces");
        axpy(1, x3, S.getPositionModifiable());
        S.getMeanCurv();
        COUT("Number of surfaces: "<<S.getNumberOfSurfaces());

        //Decreasing the number of vesicles
        S.getPositionModifiable().resize(nv2);
        ASSERT(S.getNumberOfSurfaces()==nv2,"Expeced number of surfaces");
        axpy(1, x2, S.getPositionModifiable());
        axpy(-1, S.getMeanCurv(), H2, H2);
        ASSERT(MaxAbs(H2)<1e-14,"correct computation of curvature");
        COUT(emph<<"Number of surfaces: "<<S.getNumberOfSurfaces()<<", error in curvature: "<<MaxAbs(H2)<<emph);

        //Decreasing the number of vesicles
        S.getPositionModifiable().resize(nv1);
        ASSERT(S.getNumberOfSurfaces()==nv1,"Expeced number of surfaces");
        H2.resize(nv1);
        axpy(-1, S.getMeanCurv(), H1, H2);
        axpy(1, x1, S.getPositionModifiable());
        axpy(-1, S.getMeanCurv(), H1, H1);
        ASSERT(MaxAbs(H1)<1e-14,"correct computation of curvature");
        COUT(emph<<"Number of surfaces: "<<S.getNumberOfSurfaces()
            <<", error in curvature (before update): "<<MaxAbs(H2)
            <<", error in curvature (after update): "<<MaxAbs(H1)<<emph);
    }

    //Checking the repartitioning code replicating mpi with threads
    {
        int nThreads(3);
        Repartition<real> R(&MyRepart, NULL, nThreads);
        int p(4);
        int inc(1000);

#pragma omp parallel num_threads(nThreads)
        {
            int thread(omp_get_thread_num());
            int nv = (thread + 1) * inc;
            Vec_t x(nv, p);
            Sca_t tension(nv, p);

            for(int ii=0; ii<x.size();++ii)
                *(x.begin() + ii) = inc * thread * (thread + 1) * x.getSubLength()/2 + ii ;

            for(int ii=0; ii<tension.size();++ii)
                *(tension.begin() + ii) = inc * thread * (thread + 1) * tension.getSubLength()/2 + ii ;

#pragma omp critical
            {
                COUTDEBUG("Number of surfaces on thread "<<omp_get_thread_num()<<": "<<x.getNumSubs());
            }

            R(x, tension);

#pragma omp barrier
            real xHead = *(x.begin());
            real tHead = *(tension.begin());
            real xTail = *(x.end()-1);
            real tTail = *(tension.end()-1);

            for(int ii=0; ii<x.size();++ii)
                *(x.begin()+ii) -= xHead + ii;

            for(int ii=0; ii<tension.size();++ii)
                *(tension.begin()+ii) -= tHead + ii;

            ASSERT(MaxAbs(x)<1e-14,"correct repartitioning");
            ASSERT(MaxAbs(tension)<1e-14,"correct repartitioning");
        }
    }
    COUTDEBUG(emph<<"*** Repartitioning code passed ***"<<emph);
    return 0;
}
