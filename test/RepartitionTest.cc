#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "DataIO.h"
#include "Parameters.h"
#include "Surface.h"
#include "Repartition.h"

typedef float real;
#define DT CPU 

extern const Device<DT> the_device(0);

#ifndef Doxygen_skip

void MyRepart(size_t nv, size_t stride, const real* x, 
    const real* tension, size_t* nvr, 
    real** xr, real** tensionr, void*)
{
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
    typedef Scalars<real, DT, the_device> Sca_t;
    typedef Vectors<real, DT, the_device> Vec_t;
    typedef OperatorsMats<Sca_t> Mats_t;

    {
        int const p(12);
        int const nv1(1), nv2(4), nv3(4000);
    
        //Reading operators from file
        Parameters<real> sim_par;
        bool readFromFile = true;
        Mats_t mats(readFromFile, sim_par);
        
        //Creating objects
        Vec_t x1(nv1, p); x1.fillRand();
        Vec_t x2(nv2, p); x2.fillRand();
        Vec_t x3(nv3, p); x3.fillRand();

        Surface<Sca_t, Vec_t> S(x1, mats);
        Sca_t H1(nv1, p);
        axpy(1, S.getMeanCurv(), H1);
        COUT("  Number of surfaces: "<<S.getNumberOfSurfaces()<<endl);
   
        //Increasing the number of vesicles
        S.getPositionModifiable().resize(nv2);
        axpy(1, x2, S.getPositionModifiable());
        Sca_t H2(nv2, p);
        axpy(1, S.getMeanCurv(), H2);
        COUT("  Number of surfaces: "<<S.getNumberOfSurfaces()<<endl);

        //Increasing the number of vesicles
        S.getPositionModifiable().resize(nv3);
        axpy(1, x3, S.getPositionModifiable());
        S.getMeanCurv();
        COUT("  Number of surfaces: "<<S.getNumberOfSurfaces()<<endl);

        //Decreasing the number of vesicles
        S.getPositionModifiable().resize(nv2);
        axpy(1, x2, S.getPositionModifiable());
        axpy(-1, S.getMeanCurv(), H2, H2);
        COUT("  Number of surfaces: "<<S.getNumberOfSurfaces()<<endl
            <<"   Error in curvature: "<<MaxAbs(H2)<<endl);

        //Decreasing the number of vesicles
        S.getPositionModifiable().resize(nv1);
        H2.resize(nv1);
        axpy(-1, S.getMeanCurv(), H1, H2);
        axpy(1, x1, S.getPositionModifiable());
        axpy(-1, S.getMeanCurv(), H1, H1);
        COUT("  Number of surfaces: "<<S.getNumberOfSurfaces()<<endl
            <<"   Error in curvature (before update): "<<MaxAbs(H2)<<endl
            <<"   Error in curvature (after update): "<<MaxAbs(H1)<<endl<<endl);    
    }


    //Checking the repartitioning code
    {
        int nThreads(3);
        Repartition<real> R(&MyRepart);
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
                COUT("  Number of surfaces on thread "<<omp_get_thread_num()
                    <<": "<<x.getNumSubs()<<endl);
                
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
                        
#pragma omp critical
            {
                COUT("  After repartitioning\n"
                    <<"  Number of surfaces on thread "<<omp_get_thread_num()
                    <<": "<<x.getNumSubs()<<'\n'
                    <<"  Head :"<<xHead<<",\t"<<tHead<<"\n"
                    <<"  Tail :"<<xTail<<",\t" <<tTail<<"\n"
                    <<"  Deviation for ordered :"<<MaxAbs(x)<<", "<<MaxAbs(tension)<<endl);
            }
        }
    }
    return 0;
}
