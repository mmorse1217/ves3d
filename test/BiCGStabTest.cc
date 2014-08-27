#include <unistd.h>  //for sleep()
#include <sstream>

#include "Logger.h"
#include "Error.h"
#include "Device.h"

#include "Scalars.h"
// #include "Vectors.h"
#include "HelperFuns.h"
#include "BiCGStab.h"
// #include <stdlib.h>


typedef float real;

using namespace std;

#ifndef Doxygen_skip

template<typename Container>
class MatVec
{
  public:
    typedef typename Container::value_type T;
    Container diag;

    MatVec(int nf, int p)
    {
        diag.resize(nf, p);
        size_t size = diag.size();

        T* buffer = new T[size];
        for(int ii=0; ii<size;++ii)
            buffer[ii] = (ii + 1) * (ii + 1);

        diag.getDevice().Memcpy(diag.begin(), buffer,
            size * sizeof(T), MemcpyHostToDevice);

        delete[] buffer;
    }

    void operator()(Container &x, Container &ax) const
    {
        xy(diag, x, ax);
    }
};
#endif //Doxygen_skip

extern const Device<CPU> the_cpu_dev(0);

#ifdef GPU_ACTIVE
extern const Device<GPU> the_gpu_dev(0);
#endif //GPU_ACTIVE

int main(int argc, char **argv)
{
    COUT("\n ==============================\n"
        <<"  BiCGStab Test:"
        <<"\n ==============================\n");
    sleep(.5);


    typedef Scalars<real,CPU, the_cpu_dev> ScaCPU_t;

    int p = 12;
    int nfuns(1);
    ScaCPU_t x_ref(nfuns,p), b_ref(nfuns,p);
    MatVec<ScaCPU_t> Ax(nfuns, p);
    const int iter = 100;
    const int restart = 10;
    const real tol = 1e-6;

    for(int ii(0);ii<x_ref.size();++ii)
        *(x_ref.begin() + ii) = drand48();
    Ax(x_ref, b_ref);

    ostringstream cpu_o(stringstream::out);
    ostringstream gpu_o(stringstream::out);

    {
        ScaCPU_t x(nfuns, p), b(nfuns,p);
        b.getDevice().Memcpy(b.begin(), b_ref.begin(),
            b.size() * sizeof(real), MemcpyHostToDevice);

        BiCGStab<ScaCPU_t, MatVec<ScaCPU_t> > Solver;
        int iter_in(iter);
        int rsrt(restart);
        real tt(tol);

        enum BiCGSReturn ret = Solver(Ax, x, b, rsrt, iter_in, tt);

        Ax(x,b);
        axpy((real) -1.0, b_ref, b, b);

        cpu_o<<"\n  The solver returned with: ";
        cpu_o<<"\n  CPU data :";
        cpu_o<<"\n    Residual     : "<<tt;
        cpu_o<<"\n    Iter         : "<<iter_in;
        cpu_o<<"\n    True relres  : "<<sqrt(AlgebraicDot(b,b)/AlgebraicDot(b_ref,b_ref));
    }

#ifdef GPU_ACTIVE
    {
        typedef Scalars<real,GPU, the_gpu_dev> ScaGPU_t;
        MatVec<ScaGPU_t> AxGPU(nfuns, p);

        ScaGPU_t x(nfuns, p), b(nfuns,p);
        b.getDevice().Memcpy(b.begin(), b_ref.begin(),
            b.size() * sizeof(real), MemcpyHostToDevice);

        axpy(0, x, x);

        BiCGStab<ScaGPU_t, MatVec<ScaGPU_t> > Solver;
        int iter_in(iter);
        int rsrt(restart);
        real tt(tol);

        enum BiCGSReturn ret = Solver(AxGPU, x, b, rsrt, iter_in, tt);

        AxGPU(x,b);

        ScaCPU_t b_cpu(nfuns, p);
        b.getDevice().Memcpy(b_cpu.begin(), b.begin(),
            b.size() * sizeof(real), MemcpyDeviceToHost);

        axpy((real) -1.0, b_ref, b_cpu, b_cpu);

        gpu_o<<"\n\n  The solver returned with: ";
        gpu_o<<"\n  GPU data :";
        gpu_o<<"\n    Residual     : "<<tt;
        gpu_o<<"\n    Iter         : "<<iter_in;
        gpu_o<<"\n    True relres  : "<<sqrt(AlgebraicDot(b_cpu,b_cpu)/AlgebraicDot(b_ref,b_ref));
    }
#endif //GPU_ACTIVE

    cout<<cpu_o.str();
#ifdef GPU_ACTIVE
    cout<<gpu_o.str();
#endif //GPU_ACTIVE
    cout<<endl;
    sleep(.5);

   return 0;
}
